import os
import time
import glob
import math
import tempfile
import textwrap
import subprocess
import pandas as pd
from Bio import PDB
from Bio.PDB import PDBIO, Select
import warnings
import shutil
import numpy as np
from tqdm import tqdm
from gridData import Grid
from openmm.unit import *  
from openmm import LangevinIntegrator, CustomExternalForce
from openmm.app import PDBFile, Modeller, AmberPrmtopFile, PDBFile, Simulation, CutoffNonPeriodic, HBonds
from waterkit.analysis import HydrationSites, blur_map
from System_Prep.prep_sys import system_prep
from System_Prep.split_sys import run_command
from Bio.PDB import PDBParser
from scipy.spatial.distance import cdist
import mdtraj as md
import inotify.adapters
import threading

progress = 0

def extract_box_dimensions(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline()
        x = float(first_line[6:15].strip())
        y = float(first_line[15:24].strip())
        z = float(first_line[24:33].strip())
        return math.ceil(x), math.ceil(y), math.ceil(z)

def extract_gridcenter_values(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("gridcenter"):
                _, x, y, z = line.split()
                return round(float(x), 2), round(float(y), 2), round(float(z), 2)
    raise ValueError("gridcenter values not found in the file.")


def check_error(file_path, log_file):
    if not os.path.exists(file_path):
        print(f"\033[1m\033[31mError: {file_path} not found. Check the log file for details.\033[0m")
        log_to_file(f"Error: {file_path} not found.")
        return True
    return False


def replace_residue_and_chain(input_pdb, output_pdb, new_residue_name, new_chain_id):
    """
    Replace UNK residue name and X chain ID in PDB file using BioPython.
    
    Parameters:
    -----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file
    new_residue_name : str
        The desired residue name (max 3 characters)
    new_chain_id : str
        The desired chain ID (single character)
    """
    # Suppress warnings from PDBParser
    warnings.filterwarnings('ignore')
    
    # Ensure proper formatting
    new_residue_name = new_residue_name[:3].upper()
    new_chain_id = new_chain_id[0].upper() if new_chain_id else 'A'
    
    # Parse the PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb)
    
    # Store original atom element information
    atom_elements = {}
    
    # Iterate through the structure and modify
    for model in structure:
        for chain in model:
            # Change chain ID if it's X
            if chain.id == 'X':
                chain.id = new_chain_id
            
            # Change residue names if they're UNK
            for residue in chain:
                if residue.resname.strip() == 'UNK':
                    residue.resname = new_residue_name
                
                # Store original element information for each atom
                for atom in residue:
                    # Preserve the original element case
                    atom_elements[atom.serial_number] = atom.element
    
    # Write the modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    
    # Post-process the file to fix element capitalization
    fix_element_capitalization(output_pdb, atom_elements)
    
    print(f"Modified PDB saved to {output_pdb}")

def fix_element_capitalization(pdb_file, atom_elements=None):
    """
    Fix element capitalization in PDB file.
    Ensures proper case for elements like Cl, Br, etc.
    """
    # Define proper element cases
    element_cases = {
        'CL': 'Cl', 'BR': 'Br', 'MG': 'Mg', 'MN': 'Mn', 'FE': 'Fe',
        'CO': 'Co', 'NI': 'Ni', 'CU': 'Cu', 'ZN': 'Zn', 'AS': 'As',
        'SE': 'Se', 'MO': 'Mo', 'AG': 'Ag', 'CD': 'Cd', 'SN': 'Sn',
        'SB': 'Sb', 'TE': 'Te', 'BA': 'Ba', 'HG': 'Hg', 'PB': 'Pb',
        'BI': 'Bi', 'CA': 'Ca', 'NA': 'Na', 'AL': 'Al', 'SI': 'Si',
        'LI': 'Li', 'BE': 'Be', 'NE': 'Ne', 'AR': 'Ar', 'KR': 'Kr',
        'XE': 'Xe', 'RN': 'Rn', 'HE': 'He'
    }
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    corrected_lines = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # Fix atom name (columns 13-16)
            atom_name = line[12:16]
            
            # Check if it contains a known element that needs case correction
            for wrong_case, right_case in element_cases.items():
                if wrong_case in atom_name.upper():
                    # Replace in the atom name field
                    atom_name_fixed = atom_name.upper().replace(wrong_case, right_case)
                    line = line[:12] + atom_name_fixed + line[16:]
            
            # Fix element symbol (columns 77-78) if present
            if len(line) >= 78:
                element = line[76:78].strip()
                if element.upper() in element_cases:
                    element_fixed = element_cases[element.upper()].rjust(2)
                    line = line[:76] + element_fixed + line[78:]
        
        corrected_lines.append(line)
    
    with open(pdb_file, 'w') as f:
        f.writelines(corrected_lines)


def run_tleap(ligand_name, lig_resname, lig_charge, padding, total_steps):
    log_file = 'logs/2_logs_tleap_pre_waterkit'

    def log_to_file(message):
        with open(log_file, 'a') as log:
            log.write(message + '\n')
    
    log_to_file(f"Running TLEaP on: {ligand_name}.pdb")
    print(f"\033[1m\033[34mRunning TLEaP on: \033[92m{ligand_name}.pdb\033[0m")
    print("\033[1m\033[34mRunning antechamber to generate MOL2 file...\033[0m")
    print("0% | ", end="", flush=True)

    with open(log_file, 'a') as log:
        # First antechamber command
        process1 = subprocess.run(f"antechamber -i {ligand_name}_ob.pdb -fi pdb -o {ligand_name}.pdb "
                                  f"-fo pdb -ek 'maxcyc=0'",
                                  shell=True, stdout=log, stderr=log)
        
        if process1.returncode != 0:
            log_to_file("Error in first antechamber command")
            return
        
        # Second antechamber command
        process2 = subprocess.run(f"antechamber -i {ligand_name}.pdb -fi pdb -o {ligand_name}.mol2 "
                                  f"-fo mol2 -c bcc -nc {lig_charge} -m 1 -at gaff2",
                                  shell=True, stdout=log, stderr=log)
        
        if process2.returncode != 0:
            log_to_file("Error in second antechamber command")
            return
        
        start_time = time.time()
        while process2.returncode is None:
            elapsed_time = time.time() - start_time
            progress = min(int((elapsed_time / (total_steps * 0.15)) * 100), total_steps)
            print(f"\r{progress}% | {'#' * (progress // 1)}", end="", flush=True)
            time.sleep(0.25)
    
    print("\r100% | " + "#" * 100)
    print("\033[1m\033[34mAntechamber completed successfully.\033[0m")
    log_to_file("Antechamber completed successfully.")

    print("\033[1m\033[34mRunning parmchk2 to generate FRCMOD file\033[0m")
    with open(log_file, 'a') as log:
        subprocess.run(f"parmchk2 -i {ligand_name}.mol2 -f mol2 -o {ligand_name}.frcmod",
                       shell=True, stdout=log, stderr=log)
    log_to_file("parmchk2 completed successfully.")
    
    # Ensure error check function is properly defined
    if check_error(f"{ligand_name}.frcmod", log_file):
        return
    
    tleap_commands = (f"source leaprc.gaff;\n"
                      f"loadamberparams {ligand_name}.frcmod;\n"
                      f"{lig_resname} = loadmol2 {ligand_name}.mol2;\n"
                      f"check {lig_resname};\n"
                      f"setBox {lig_resname} vdw {padding};\n"  
                      f"saveamberparm {lig_resname} {ligand_name}.prmtop {ligand_name}.inpcrd;\n"
                      f"saveoff {lig_resname} {ligand_name}.lib;\n"
                      f"savepdb {lig_resname} {ligand_name}_box.pdb;\n"
                      f"quit")

    with tempfile.NamedTemporaryFile(delete=False) as tleap_file:
        tleap_file.write(tleap_commands.encode('utf-8'))
        tleap_file.close()

    with open(log_file, 'a') as log:
        subprocess.run(f"tleap -f {tleap_file.name} > logs/2_logs_tleap_pre_waterkit", shell=True, stdout=log, stderr=log)
    os.remove(tleap_file.name)

    if check_error(f"{ligand_name}.frcmod", log_file) or check_error(f"{ligand_name}.lib", log_file) or check_error(f"{ligand_name}.mol2", log_file) or check_error(f"{ligand_name}_box.pdb", log_file):
        return
    
    log_to_file(f"TLEaP run completed. Files saved: {ligand_name}.frcmod, {ligand_name}.lib, {ligand_name}.mol2, {ligand_name}_box.pdb")

    log_to_file(f"TLEaP run successfully and following files are saved:")
    log_to_file(f"{ligand_name}.frcmod, {ligand_name}.lib, {ligand_name}.mol2, {ligand_name}_box.pdb")

    print(f"\033[1m\033[34mTLEaP run successfully and following files are saved\033[0m")
    print(f"\033[1m\033[34m\033[92m{ligand_name}.frcmod {ligand_name}.lib {ligand_name}.mol2 {ligand_name}_box.pdb\033[0m\n\n")

    return f"{ligand_name}.frcmod, {ligand_name}.lib, {ligand_name}.mol2, {ligand_name}_box.pdb"


def process_pdb(input_pdb, lig_resname=None):
    log_file = 'logs/1_logs_process_pdb'

    def log_to_file(message):
        with open(log_file, 'a') as log:
            log.write(message + '\n')

    os.makedirs("logs", exist_ok=True)

    with open(input_pdb, 'r') as file:
        lines = file.readlines()

    filtered_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM', 'TER'))]

    antechamber_output_ligand_pdb = 'ligand_wH.pdb'
    with open(antechamber_output_ligand_pdb, 'r') as file:
        lines = file.readlines()

    filtered_ligand_lines = [line for line in lines if line.startswith('ATOM') and lig_resname in line]

    start_index = None
    end_index = None
    for i, line in enumerate(filtered_lines):
        if line.startswith('HETATM') and lig_resname in line:
            if start_index is None:
                start_index = i
            end_index = i

    if start_index is not None and end_index is not None:
        modified_lines = filtered_lines[:start_index] + filtered_ligand_lines + filtered_lines[end_index + 1:]
    else:
        modified_lines = filtered_lines + filtered_ligand_lines

    output_pdb = input_pdb.replace('.pdb', '_wH.pdb')
    with open(output_pdb, 'w') as file:
        file.writelines(modified_lines)
        file.write('END\n')

    if not os.path.exists(output_pdb):
        error_message = "Error: Output PDB file creation failed! Check the log file for more details."
        print(f"\033[1m\033[31m{error_message}\033[0m")
        log_to_file(error_message)
        return

    log_to_file(f"Protein PDB with added hydrogens in ligand saved as: {output_pdb}")
    print(f"\n\033[1m\033[34mProtein PDB with added hydrogens in ligand saved as: \033[92m{output_pdb}\033[0m")

    return output_pdb
    
def prepare_receptor_and_create_grid_protein(SX, SY, SZ, input_pdb, output_pdb, ligand_name, water_model, skip_h_add):

    log_file = 'logs/3_logs_prepare_receptor_and_create_grid_protein'

    def log_to_file(message):
        with open(log_file, 'a') as log:
            log.write(message + '\n')
            
    pdb_file = input_pdb if skip_h_add else output_pdb
    print(pdb_file)
    
    with open(log_file, 'a') as log:
        subprocess.run(f"wk_prepare_receptor.py -i {pdb_file} --lib {ligand_name}.lib --frcmod {ligand_name}.frcmod \
                         -o prot_prep --ignore_gaps --pdb --amber_pdbqt",
                         shell=True, stdout=log, stderr=log)
        
        subprocess.run(f"/root/../miniconda/envs/waterkit/bin/wk_create_grid_protein_file.py \
                         -r prot_prep_amber.pdbqt -l {ligand_name}.pdb \
                         -s {SX} {SY} {SZ} -w \"{water_model}\" -o prot_lig.gpf", 
                         shell=True, stdout=log, stderr=log)

        log_to_file(f"Box Dimensions: {SX}, {SY}, {SZ}")
        log_to_file(f"Receptor and grid preparation completed successfully.")
        print(f"\033[1m\033[34mReceptor and grid preparation completed successfully.\033[0m")


def check_for_errors(log_file):
    """Check the log file for errors."""
    error_keywords = ["Traceback", "KeyError", "ValueError", "Cannot find", "Error"]
    with open(log_file, 'r') as log:
        log_content = log.read()
        return any(keyword in log_content for keyword in error_keywords)

def run_waterkit(X, Y, Z, SX, SY, SZ, n_cores, total_steps):
    log_file = 'logs/4_logs_waterkit'

    def log_to_file(message):
        with open(log_file, 'a') as log:
            log.write(message + '\n')
            
    os.makedirs("traj", exist_ok=True)
    
    log_to_file(f"Running WaterKit simulation...")
    print("\n\033[1m\033[34mRunning WaterKit simulation...\033[0m")
    print("0% | ", end="", flush=True)
    
    process = subprocess.Popen(f'run_waterkit.py -i prot_prep_amber.pdbqt \
                                 -c {X} {Y} {Z} -s {SX} {SY} {SZ} -n {total_steps} -j {n_cores} -o traj',
                                 shell=True, stdout=open(log_file, 'a'), stderr=open(log_file, 'a'))
    
    while process.poll() is None:
        num_files = len([name for name in os.listdir('traj') if os.path.isfile(os.path.join('traj', name))])
        progress = min(int((num_files / total_steps) * 100), 100)
        print(f"\r{progress}% | {'#' * progress}", end="", flush=True)
        time.sleep(1)
    print("\r100% | " + "#" * 100)
    
    if check_for_errors(log_file):
        print("\033[1m\033[31mError detected in WaterKit simulation! Please check log file 4_logs_waterkit.\033[0m")
        return False  

    log_to_file(f"Making trajectories...")
    print("\n\n\033[1m\033[34mMaking trajectories...\033[0m")
    print("0% | ", end="", flush=True)
    
    process = subprocess.Popen('/root/../miniconda/envs/waterkit/bin/wk_make_trajectory.py \
                               -r prot_prep.pdb -w traj -o prot_lig_out --lib ligand_wH.lib --frcmod ligand_wH.frcmod',
                               shell=True, stdout=open(log_file, 'a'), stderr=open(log_file, 'a'))
    
    while process.poll() is None:
        time.sleep(1)
    print("\r100% | " + "#" * 100)
    
    if check_for_errors(log_file):
        print("\033[1m\033[31mError detected in trajectory creation! Please check log file 4_logs_waterkit.\033[0m")
        return False

    log_to_file(f"Minimizing trajectories...")
    print("\033[1m\033[34mMinimizing trajectories...\033[0m")
    print("0% | ", end="", flush=True)
    
    process = subprocess.Popen(f'wk_minimize_trajectory.py -p prot_lig_out_system.prmtop \
                                 -t prot_lig_out_system.nc -o prot_lig_out_system_min.nc',
                                 shell=True, stdout=open(log_file, 'a'), stderr=open(log_file, 'a'))
    
    while process.poll() is None:
        time.sleep(0.25)
    print("\r100% | " + "#" * 100)
    
    if check_for_errors(log_file):
        print("\033[1m\033[31mError detected in minimization! Please check log file 4_logs_waterkit.\033[0m")
        return False
    
    print("\033[1m\033[34mMinimizing trajectories completed.\033[0m\n\n")
    return True  


def run_gist(X, Y, Z, SX, SY, SZ, gist_res):
    log_file = 'logs/5_logs_gist'

    def log_to_file(message):
        with open(log_file, 'a') as log:
            log.write(message + '\n')

    with open(log_file, 'a') as log:
        log_to_file(f"Running GIST analysis...")
        print("\033[1m\033[34mRunning GIST analysis...\033[0m")
        grid_factor = int(1 / gist_res)
        gist_input_content = f"""parm prot_lig_out_system.prmtop
trajin prot_lig_out_system_min.nc
gist gridspacn {gist_res} gridcntr {X} {Y} {Z} griddim {grid_factor * SX} {grid_factor * SY} {grid_factor * SZ}
go
quit
    """
    
        with open('gist.inp', 'w') as gist_file:
            gist_file.write(gist_input_content)
    
        process = subprocess.Popen("cpptraj -i gist.inp", shell=True, stdout=log, stderr=log)
        time.sleep(2)

        log_to_file(f"GIST analysis completed successfully")
        print("\033[1m\033[34mGIST analysis completed successfully and following files saved successfully.\033[0m")
        
        gist_files = glob.glob('gist-*.dx')
        relevant_gist_files = [f for f in gist_files if any(f.startswith(prefix) for prefix in ['gist-gO', 'gist-Esw', 'gist-Eww', 'gist-dTStrans', 'gist-dTSorient'])]
        
#        if not relevant_gist_files:
#            print("\033[1m\033[31mError: No relevant GIST files generated! Please check log file 5_logs_gist. \033[0m")
#            return  
        
        print("\033[1m\033[34m\033[92m" + " ".join(relevant_gist_files) + "\033[0m")

        return True


def parse_pdb_line(line):
    """Parse a PDB ATOM/HETATM line with detailed debugging"""
    if not line.startswith(('ATOM', 'HETATM')):
        return None
    
    try:
        # Handle chain ID that might be blank
        chain_id = line[21:22]
        if chain_id == ' ':
            chain_id = ''  # Convert space to empty string
        
        return {
            'record': line[0:6].strip(),
            'atom_num': int(line[6:11]),
            'atom_name': line[12:16].strip(),
            'resname': line[17:20].strip(),
            'chain': chain_id,
            'resnum': int(line[22:26]),
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54])
        }
    except (ValueError, IndexError) as e:
        print(f"Error parsing line: {line.strip()}")
        print(f"Error: {e}")
        return None

def find_closest_water(crystal_pdb, target_pdb, resid_crystal):
    
    # Load PDB structures
    crystal = md.load_pdb(crystal_pdb)
    prot_opmm = md.load_pdb(target_pdb)

    # Parse residue ID (handle both '1326' and 'A1326' formats)
    def parse_resid(resid_str):
        resid_str = str(resid_str)
        if resid_str.isdigit():
            return None, int(resid_str)
        else:
            chain = resid_str[0]
            resnum = int(resid_str[1:])
            return chain, resnum

    chain, resnum = parse_resid(resid_crystal)

    # Use MDTraj selection language to find the water oxygen directly
    if chain is None:
        # No chain specified - search by resnum only
        selection_string = f"resname HOH and name O and resSeq {resnum}"
    else:
        # Chain specified - include chain in selection
        # Convert chain letter to chain index (A=0, B=1, etc.)
        chain_index = ord(chain.upper()) - ord('A')
        selection_string = f"resname HOH and name O and resSeq {resnum} and chainid {chain_index}"
    
    # Get the oxygen atom index directly
    oxygen_atoms = crystal.topology.select(selection_string)
    
    if len(oxygen_atoms) == 0:
        print(f"Residue {resid_crystal} not found in {crystal_pdb}!")
        # Show available water residues for debugging
        water_residues = crystal.topology.select("resname HOH")
        if len(water_residues) > 0:
            print("Available water residues:")
            seen_residues = set()
            for idx in water_residues:
                atom = crystal.topology.atom(idx)
                res_info = f"{atom.residue.chain.index}:{atom.residue.resSeq}"
                if res_info not in seen_residues:
                    seen_residues.add(res_info)
                    print(f"  Chain {chr(atom.residue.chain.index + ord('A'))}, ResID {atom.residue.resSeq}")
        return None
    
    oxygen_atom_index = oxygen_atoms[0]  # Take the first match if multiple

    # Step 3: Get coordinates of the selected water oxygen
    water_crystal = crystal.xyz[0][oxygen_atom_index] * 10  # Convert nm → Å

    # Step 4: Get all water oxygens in target_pdb
    water_indices_opmm = prot_opmm.topology.select("resname HOH and name O")

    if len(water_indices_opmm) == 0:
        print(f"No water molecules found in {target_pdb}!")
        return None

    # Step 5: Extract coordinates of all water oxygens in target_pdb
    water_coords = prot_opmm.xyz[0][water_indices_opmm] * 10  # Convert nm → Å

    # Step 6: Compute distances between crystal water oxygen and target_pdb waters
    distances = np.linalg.norm(water_coords - water_crystal, axis=1)

    # Step 7: Find the closest water molecule
    closest_index = water_indices_opmm[np.argmin(distances)]
    closest_resid = prot_opmm.topology.atom(closest_index).residue.resSeq
    
    # Get the chain of the closest water if available
    closest_chain = prot_opmm.topology.atom(closest_index).residue.chain
    closest_chain_id = chr(closest_chain.index + ord('A'))

    print(f"Closest water in {target_pdb} is chain {closest_chain_id} resid {closest_resid} at distance {np.min(distances):.3f} Å")
    
    # Return with chain if it was specified in input
    if chain is not None:
        new_wat_site = f"{closest_chain_id}{closest_resid}"
    else:
        new_wat_site = str(closest_resid)
    
    return new_wat_site

def find_nearby_waters(lig_resname, wat_chain_id="A", wat_sel_dist=5.0):
    ligand_coords = []
    water_coords = []
    water_resids = []
    
    ligand_count = 0
    water_count = 0
    total_lines = 0
    
    print(f"Looking for ligand: {lig_resname}")
    print(f"Looking for waters in chain: '{wat_chain_id}'")
    print(f"Distance threshold: {wat_sel_dist} Å")
    print("-" * 50)
    
    with open('prot_lig_out_system_wk.pdb', 'r') as f:
        for line_num, line in enumerate(f, 1):
            total_lines += 1
            atom_data = parse_pdb_line(line)
            if not atom_data:
                continue
            
            # Debug: Print first few atoms to see what we're getting
            if total_lines <= 10:
                print(f"Line {line_num}: {atom_data}")
            
            # Collect ligand atoms
            if atom_data['resname'] == lig_resname:
                if not atom_data['atom_name'].startswith('H'):
                    ligand_coords.append([atom_data['x'], atom_data['y'], atom_data['z']])
                    ligand_count += 1
                    if ligand_count <= 5:  # Print first few ligand atoms
                        print(f"LIGAND atom {ligand_count}: {atom_data}")
            
            # Collect water atoms - let's be very explicit about the conditions
            is_water_resname = atom_data['resname'] in ['WAT', 'HOH', 'SOL']
            is_correct_chain = atom_data['chain'] == wat_chain_id
            
            if is_water_resname and is_correct_chain:
                water_coords.append([atom_data['x'], atom_data['y'], atom_data['z']])
                water_resids.append(f"{atom_data['chain']}{atom_data['resnum']}")
                water_count += 1
                if water_count <= 5:  # Print first few water atoms
                    print(f"WATER atom {water_count}: {atom_data}")
            elif is_water_resname:
                # Debug: show waters in wrong chain
                if water_count == 0:  # Only show this once
                    print(f"Found water in wrong chain: chain='{atom_data['chain']}', expected='{wat_chain_id}'")
    
    print("-" * 50)
    print(f"Total lines processed: {total_lines}")
    print(f"Found {ligand_count} ligand atoms")
    print(f"Found {water_count} water atoms")
    print(f"Found {len(set(water_resids))} unique water residues")
    
    if ligand_count == 0:
        raise ValueError(f"No ligand atoms found with residue name {lig_resname}")
    
    if water_count == 0:
        # Let's see what chains we actually have for waters
        print("\nDebugging: Let's see all water residues and their chains...")
        with open('prot_lig_out_system_wk.pdb', 'r') as f:
            water_chains_found = set()
            for line in f:
                atom_data = parse_pdb_line(line)
                if atom_data and atom_data['resname'] in ['WAT', 'HOH', 'SOL']:
                    water_chains_found.add(f"'{atom_data['chain']}'")
                    if len(water_chains_found) <= 3:  # Show first few
                        print(f"Water found in chain: '{atom_data['chain']}'")
            print(f"All water chains found: {sorted(water_chains_found)}")
        raise ValueError(f"No water molecules found in chain '{wat_chain_id}'")
    
    # Calculate distances
    ligand_coords = np.array(ligand_coords)
    water_coords = np.array(water_coords)
    
    print(f"\nLigand coords shape: {ligand_coords.shape}")
    print(f"Water coords shape: {water_coords.shape}")
    
    distances = cdist(ligand_coords, water_coords)
    print(f"Distance matrix shape: {distances.shape}")
    print(f"Min distance: {distances.min():.2f} Å")
    print(f"Max distance: {distances.max():.2f} Å")
    
    within_threshold = (distances < wat_sel_dist).any(axis=0)
    nearby_water_resids = list(set(np.array(water_resids)[within_threshold]))
    
    print(f"Waters within {wat_sel_dist} Å: {len(nearby_water_resids)}")
    print(f"Nearby water resids within {wat_sel_dist} Å of ligand ({lig_resname}):")
    print(sorted(nearby_water_resids))
    
    if nearby_water_resids:
        return sorted(nearby_water_resids)
    else:
        raise ValueError(f"No water residue found within {wat_sel_dist} Å of the ligand.")


def generate_water_site_pdb(wat_sites, input_pdb):
    water_site_lines = {site: [] for site in wat_sites}

    with open(input_pdb, "r") as input_file:
        for line in input_file:
            if (line.startswith("HETATM") or line.startswith("ATOM")) and any(w in line for w in ["HOH", "WAT", "SOL"]):
                chain_id = line[21]
                residue_number = line[22:26].strip()
                residue_id = f"{chain_id}{residue_number}"
                if residue_id in water_site_lines:
                    water_site_lines[residue_id].append(line)

    for wat_site, lines in water_site_lines.items():
        with open(f"wat_site_{wat_site}.pdb", "w") as output_file:
            output_file.writelines(lines)

def create_summary_file(log_file_path, output_file_path='../wat_FE_summary.dat'):
    """
    Reads a gistpp.log file, extracts free energy values for water residues,
    and writes them to a wat_FE_summary.dat file.
    
    Parameters:
    log_file_path (str): Path to the input gistpp.log file.
    output_file_path (str): Path to the output wat_FE_summary.dat file.
    """
    try:
        with open(log_file_path, 'r') as file:
            lines = file.readlines()

        data = []
        current_site = None
        free_energy = None
        tds = None
        etot = None

        for line in lines:
            if "Processing wat_site:" in line:
                current_site = line.split(":")[1].strip()
                # Reset values for new site
                free_energy = None
                tds = None
                etot = None
            elif "sum of: wat_dA_" in line and "is:" in line:
                # Extract value after "is:"
                value_str = line.split("is:")[1].strip()
                free_energy = float(value_str)
            elif "sum of: wat_dTStot_" in line and "is:" in line:
                value_str = line.split("is:")[1].strip()
                tds = float(value_str)
            elif "sum of: wat_Etot_" in line and "is:" in line:
                value_str = line.split("is:")[1].strip()
                etot = float(value_str)
            elif "Completed processing wat_site:" in line:
                if current_site is not None and all(v is not None for v in [free_energy, tds, etot]):
                    data.append((current_site, free_energy, tds, etot))
                current_site = None

        # Sort data by water site name for consistent output
        data.sort(key=lambda x: x[0])

        with open(output_file_path, 'w') as file:
            file.write("Water_Site  FreeEnergy_in_kcal/mol  TdS_in_kcal/mol  Etot_in_kcal/mol\n")
            for site, free_energy, tds, etot in data:
                file.write(f"{site:<10}  {free_energy:>20.6f}  {tds:>16.6f}  {etot:>17.6f}\n")
                    
        print(f"Summary file created successfully at: {output_file_path}")
        print(f"Processed {len(data)} water sites")
        
    except FileNotFoundError:
        print(f"Error: The file {log_file_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def add_water_sites(filename, old_water_sites_list):
    """
    Add a 'Crystal_water_site' column as the first column in the file.
    Each row gets assigned a value from old_water_sites_list sequentially.
    
    Parameters:
    - filename: Input file name
    - old_water_sites_list: List of crystal water sites to assign to rows
    """
    
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Process the file
    output_lines = []
    
    # Process header line - add Crystal_water_site as first column
    header_parts = lines[0].strip().split()
    new_header = "Crystal_water_site\t" + "\t".join(header_parts) + "\n"
    output_lines.append(new_header)
    
    # Process data lines
    data_line_index = 0
    for line in lines[1:]:
        line = line.strip()
        if line:  # Skip empty lines
            # Get the crystal water site for this row (if available)
            if data_line_index < len(old_water_sites_list):
                crystal_site = old_water_sites_list[data_line_index]
            else:
                crystal_site = ""  # Empty if we run out of sites
            
            # Split the line and reconstruct with new first column
            parts = line.split()
            new_line = crystal_site + "\t" + "\t".join(parts) + "\n"
            output_lines.append(new_line)
            data_line_index += 1
    
    # Overwrite the original file
    with open(filename, 'w') as f:
        f.writelines(output_lines)
    
    print(f"File '{filename}' has been updated with Crystal_water_site column")

def run_gistpp(lig_resname, wat_chain_id="A", wat_sel_dist=5.0, wat_sites=None):
    if wat_sites:
        generate_water_site_pdb(wat_sites,'prot_lig_out_system_wk.pdb')
    else:
        wat_sites = find_nearby_waters(lig_resname, wat_chain_id, wat_sel_dist)
        if wat_sites is not None:
            generate_water_site_pdb(wat_sites,'prot_lig_out_system_wk.pdb')
    command = f"bash gistpp.sh {' '.join(wat_sites)}"
    run_command(command, "gistpp.log")

    create_summary_file('gistpp.log', 'wat_FE_summary.dat')



def hydra_sites(gist_res):
    log_file = 'logs/6_logs_hydra_sites'
    
    def log_to_file(message):
        os.makedirs("logs", exist_ok=True)  
        with open(log_file, 'a') as log:
            log.write(message + '\n')

    with open(log_file, 'a') as log:
        required_files = {"gist-gO.dx", "gist-Esw-dens.dx", "gist-Eww-dens.dx", "gist-dTStrans-dens.dx", "gist-dTSorient-dens.dx"}
        
        log_to_file(f"Waiting for GIST files to be generated...")
        print("\033[1m\033[34mWaiting for GIST files to be generated...\033[0m")
        
        # First, check if files already exist
        found_files = {f for f in required_files if os.path.exists(f)}
        if found_files:
            print(f"Already found: {', '.join(found_files)}")
        
        # Only watch if we need to wait for files
        if found_files != required_files:
            i = inotify.adapters.Inotify()
            i.add_watch('.')
            
            for event in i.event_gen(yield_nones=False):
                (_, type_names, path, filename) = event
                # Watch for multiple events that indicate file creation/completion
                if filename in required_files and any(event_type in type_names for event_type in ['IN_CREATE', 'IN_CLOSE_WRITE', 'IN_MOVED_TO']):
                    if os.path.exists(filename):  # Double-check the file actually exists
                        found_files.add(filename)
                        print(f"Found: {filename}")
                        if found_files == required_files:
                            break
        
        log_to_file(f"All required files found. Starting hydration site identification...")
        print("\033[1m\033[34mAll required files found. Starting hydration site identification...\033[0m")
        
        gO = Grid("gist-gO.dx")
        esw = Grid('gist-Esw-dens.dx')
        eww = Grid('gist-Eww-dens.dx')
        tst = Grid('gist-dTStrans-dens.dx')
        tso = Grid('gist-dTSorient-dens.dx')
        dg = (esw + 2 * eww) - (tst + tso)
    
        hs = HydrationSites(gridsize=gist_res, water_radius=1.4, min_water_distance=2.5, min_density=0.1)        
        hydration_sites = hs.find(gO)
        dg_energy = hs.hydration_sites_energy(dg, water_radius=1.4)
    
        print("0% | ", end="", flush=True)
        total_steps = 100
        update_interval = 1
        start_time = time.time()
    
        while time.time() - start_time < total_steps * update_interval:
            elapsed_time = time.time() - start_time
            progress = min(int((elapsed_time / (total_steps * update_interval)) * 100), 100)
            print(f"\r{progress}% | {'#' * progress}", end="", flush=True)
            time.sleep(update_interval)
    
        hs.export_to_pdb("water_sites_dG_smoothed.pdb", hydration_sites, dg_energy)
        dg_energy_df = pd.DataFrame(dg_energy, columns=["Free Energy (~TG)"])
        dg_energy_df.to_csv("water_sites_dG_energy.txt", index=False, header=False, float_format="%.3f")
        map_smooth = blur_map(dg, radius=1.4)
        map_smooth.export("gist_dG_dens_smoothed.dx")
        
        print("\r100% | " + "#" * 100)
        log_to_file(f"Identification of hydration sites completed.")
        print("\033[1m\033[34mIdentification of hydration sites completed.\033[0m\n")

        return True


def process_wk_pdb():
    with open("prot_lig_out_system.pdb", 'r') as file:
        lines = file.readlines()
    filtered_lines = [line for line in lines if not (line.startswith('HETATM') and 'HOH' in line) and 'END' not in line]

    with open('water_sites_dG_smoothed.pdb', 'r') as file:
        water_lines = file.readlines()
    water_lines = [line.replace('HETATM', 'ATOM  ') for line in water_lines]
    
    modified_lines = filtered_lines + water_lines

    output_pdb = "prot_lig_out_system.pdb".replace('.pdb', '_wk.pdb')
    with open(output_pdb, 'w') as file:
        file.writelines(modified_lines)
        file.write('END\n')

    print(f"\033[1m\033[34mprot_lig_out_system.pdb and water_sites_dG_smoothed.pdb files are merged and saved as:\033[0m \033[1m\033[34m\033[92m{output_pdb}\033[0m")


def add_hydrogens():
    pdb = PDBFile('water_sites_dG_smoothed.pdb')
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens()

    with open('water_sites_wH.pdb', 'w') as file:
        PDBFile.writeFile(modeller.topology, modeller.positions, file)

    with open("prot_lig_out_system.pdb", 'r') as file:
        lines = file.readlines()

    filtered_lines = []
    for line in lines:
        if line.startswith('HETATM') and 'HOH' in line:
            continue
        if line.startswith('END'):
            continue
        filtered_lines.append(line)

    with open('water_sites_wH.pdb', 'r') as file:
        h_lines = file.readlines()

    h_lines = [line for line in h_lines if line.startswith('HETATM')]

    modified_lines = filtered_lines + h_lines

    output_pdb = "prot_lig_wat_modH.pdb"
    with open(output_pdb, 'w') as file:
        file.writelines(modified_lines)
        file.write('END\n')

    print(f"\033[1m\033[34mFinal pdb file is saved as:\033[0m \033[1m\033[34m\033[92m{output_pdb}\033[0m")

def tleap_post_wkit_comp(lig_resname, water_model):

    tleap_commands = (f"source leaprc.protein.ff14SB;\n"
                      f"source leaprc.water.{water_model};\n"
                      f"source leaprc.gaff2;\n"
                      f"loadamberparams ligand_wH.frcmod;\n"
                      f"{lig_resname} = loadmol2 ligand_wH.mol2;\n"
                      f"check {lig_resname};\n"
                      f"comp = loadpdb prot_lig_wat_modH.pdb;\n"
                      f"saveamberparm comp prot_lig_wat_modH.prmtop prot_lig_wat_modH.inpcrd;\n"
                      f"savepdb comp prot_lig_wat_fin_tleap.pdb;\n" 
                      f"quit")
    
    with tempfile.NamedTemporaryFile(delete=False) as tleap_file:
        tleap_file.write(tleap_commands.encode('utf-8'))
        tleap_file.close()  
    
    subprocess.run(f"tleap -f {tleap_file.name} > logs/7_logs_tleap_post_waterkit", shell=True)
    os.remove(tleap_file.name)

def openmm_min():
    prmtop_file = "prot_lig_wat_modH.prmtop"
    pdb_file = "prot_lig_wat_modH.pdb"
    output_pdb = "prot_opmm_min.pdb"

    prmtop = AmberPrmtopFile(prmtop_file)
    pdb = PDBFile(pdb_file)

    system = prmtop.createSystem(nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=9*angstroms, constraints=HBonds)

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    print("\n\033[1m\033[34mRunning OpenMM...\033[0m")
    simulation.minimizeEnergy(tolerance=1*kilojoule/(nanometer*mole), maxIterations=0)
    state = simulation.context.getState(getEnergy=True)
    min_energy = state.getPotentialEnergy()

    energy_log = "min_energy.log"
    with open(energy_log, "w") as file:
        file.write(f"Minimized energy: {min_energy}\n")

    with open(output_pdb, 'w') as output:
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), output)
    print("\033[1m\033[34mOpenMM calculation completed.\033[0m")
        
def process_pdb_and_run_tleap(input_pdb, prot_chain_id=None, lig_chain_id=None, lig_resid=None, lig_resname=None, 
                                lig_charge=None, padding=None, water_model=None, wk_frames=None, n_cores=None, 
                                gist_res=None, skip_h_add=False, wat_chain_id="A", wat_sel_dist=5.0, wat_sites=None):
    start_time = time.time()
    total_steps = wk_frames
    os.makedirs("logs", exist_ok=True)
    system_prep(input_pdb, prot_chain_id, lig_chain_id, lig_resid, lig_resname)
    input_pdb_crys = input_pdb
    replace_residue_and_chain("prep_min_lig.pdb", "ligand_wH_ob.pdb", lig_resname, lig_chain_id)
    if not run_tleap("ligand_wH", lig_resname, lig_charge, padding, total_steps):
        print("\033[1m\033[31mTleap failed. Please check log file 2_logs_run_tleap.\033[0m")
        return

    SX, SY, SZ = extract_box_dimensions("ligand_wH_box.pdb")
    print(f"Box Dimensions: {SX}, {SY}, {SZ}")

    input_pdb = "prep_min_prot_wat.pdb"
    output_pdb = process_pdb(input_pdb, lig_resname)
    
    if not output_pdb:
        print("\033[1m\033[31mProcessing pdb failed. Please check log file 1_logs_process_pdb.\033[0m")
        return

    prepare_receptor_and_create_grid_protein(SX, SY, SZ, input_pdb, output_pdb, "ligand_wH", water_model, skip_h_add)

    X, Y, Z = extract_gridcenter_values("prot_lig.gpf")
    print(f"Grid Centers (Centroids): {X}, {Y}, {Z}")

    if not run_waterkit(X, Y, Z, SX, SY, SZ, n_cores, total_steps):
        return
    if not run_gist(X, Y, Z, SX, SY, SZ, gist_res):
        return
    if not hydra_sites(gist_res):
        return
    process_wk_pdb()
    if len(input_pdb_crys) == 4 and not os.path.isfile(input_pdb_crys):
        crystal_pdb = f"{input_pdb_crys.upper()}.pdb"
    else:
        crystal_pdb = input_pdb_crys
    if wat_sites:
        new_wat_sites = []
        for wat_site in wat_sites:
            new_wat_site = find_closest_water(crystal_pdb,'prot_lig_out_system_wk.pdb',wat_site[1:])
            new_wat_sites.append(f"{wat_site[0]}{new_wat_site}")
        print(new_wat_sites)
        run_gistpp(lig_resname, wat_chain_id=wat_chain_id, wat_sel_dist=wat_sel_dist, wat_sites=new_wat_sites)
        add_water_sites('wat_FE_summary.dat',wat_sites)
    else:
        run_gistpp(lig_resname, wat_chain_id=wat_chain_id, wat_sel_dist=wat_sel_dist, wat_sites=None)
    add_hydrogens()
    tleap_post_wkit_comp(lig_resname, water_model)
    openmm_min()

    for pattern in ["ANTECHAMBER*", "sqm*", "*.inpcrd", "*.log", "*.prmtop", "*.INF", "gist-dTSsix-dens.dx", 
                    "gist-dipole*", "gist-gH.dx", "*norm.dx", "gist-output.dat", "gist_dG_dens_smoothed.dx"]:
        for file_path in glob.glob(pattern):
            os.remove(file_path)
    
    end_time = time.time()  
    total_time = (end_time - start_time) / 60 
    with open("logs/total_execution_time.txt", "w") as file:
        file.write(f"Total Execution Time: {total_time:.2f} minutes\n")
    print(f"\n\033[1m\033[34mTotal Execution Time: \033[0m \033[1m\033[34m\033[92m{total_time:.2f} minutes\033[0m\n")


