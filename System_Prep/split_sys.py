import os
import sys
import shutil
import subprocess
import requests
from Bio.PDB import PDBParser, PDBIO, Select
import re
from typing import List, Tuple, Dict
from bs4 import BeautifulSoup
import math
from rdkit import Chem

def check_hydrogens_pdb(pdb_file: str):
    # Load ligand from PDB with explicit hydrogens kept
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not parse {pdb_file}")

    total_H = 0
    polar_H = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # hydrogen
            total_H += 1
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() in (7, 8, 16):  # N, O, S
                polar_H += 1

    # classify
    if total_H == 0:
        status = "no explicit hydrogens"
    elif total_H == polar_H:
        status = "only polar hydrogens"
    else:
        status = "all hydrogens present"

    return {
        "status": status,
        "total_H": total_H,
        "polar_H": polar_H
    }


def download_sdf_files(pdb_id, ligand_resname):
    # 1. Get the structure page
    url = f'https://www.rcsb.org/structure/{pdb_id}'
    r = requests.get(url)
    soup = BeautifulSoup(r.text, 'html.parser')

    # 2. Find the Ligands Table
    ligand_table = soup.find('table', {'id': 'LigandsMainTable'})
    if ligand_table is None:
        raise Exception("Could not find Ligands table on page!")

    # 3. Find the row(s) for specified ligand
    ligand_row = ligand_table.find('tr', {'id': f'ligand_row_{ligand_resname}'})
    if ligand_row is None:
        raise Exception(f"Ligand {ligand_resname} not found in ligand table for PDB {pdb_id}.")

    # 4. Look for all SDF download URLs in the dropdown menu
    sdf_urls = []
    for a in ligand_row.find_all('a', href=True):
        href = a['href']
        if 'encoding=sdf' in href and href.endswith('.sdf'):
            # Full absolute URL
            if not href.startswith('http'):
                href = 'https://models.rcsb.org' + href
            sdf_urls.append(href)

    print(f"Found SDF URLs: {sdf_urls}")

    if not sdf_urls:
        print("No SDF file links found.")
        return

    filenames = []
    # 5. Download the SDF files
    for url in sdf_urls:
        filename =  url.split('=')[-1]
        filenames.append(filename)
        print(f"Downloading {url} ...")
        resp = requests.get(url)
        if resp.ok:
            with open(filename, 'wb') as f:
                f.write(resp.content)
            print(f"Saved: {filename}")
        else:
            print(f"Failed to download {url}")
    return filenames


def parse_pdb_atoms(pdb_file):
    atoms = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM  "):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = line[76:78].strip()
                atoms.append((x, y, z, element))
    return atoms

def parse_sdf_atoms(sdf_file):
    atoms = []
    with open(sdf_file) as f:
        lines = f.readlines()
        if len(lines) < 4:
            return []
        counts_line = lines[3]
        atom_count = int(counts_line[0:3].strip())
        atom_lines = lines[4:4 + atom_count]
        for line in atom_lines:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
            element = line[31:34].strip()
            atoms.append((x, y, z, element))
    return atoms

def atom_is_matched(s_atom, pdb_atoms, tolerance=0.15):
    xs, ys, zs, es = s_atom
    for (xp, yp, zp, ep) in pdb_atoms:
        dist = math.sqrt((xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2)
        if dist < tolerance and ep == es:
            return True
    return False

def significant_overlap(pdb_atoms, sdf_atoms, min_overlap=0.8):
    if not sdf_atoms or not pdb_atoms:
        return False
    matched = 0
    for s_atom in sdf_atoms:
        if atom_is_matched(s_atom, pdb_atoms):
            matched += 1
    overlap_ratio = matched / len(sdf_atoms)
    return overlap_ratio >= min_overlap

def get_matching_sdf_filename(pdb_filename, sdf_filename, min_overlap=0.8, tolerance=0.15):
    pdb_atoms = parse_pdb_atoms(pdb_filename)
    sdf_atoms = parse_sdf_atoms(sdf_filename)
    if significant_overlap(pdb_atoms, sdf_atoms, min_overlap=min_overlap):
        return sdf_filename
    return None

def calculate_charge(atom_symbol: str, bond_count: int) -> int:
    """
    Calculate formal charge based on:
    Formal charge = Valence electrons - Non-bonding electrons - (Bonding electrons / 2)
    
    Since we don't have lone pair info, we use:
    Formal charge = Normal valence - Bond order
    """
    
    # Normal valence (group valence) for neutral atoms
    normal_valence = {
        'H': 1,
        'C': 4,
        'N': 3,
        'O': 2,
        'P': 3,
        'S': 2,
        'F': 1,
        'Cl': 1,
        'Br': 1,
        'I': 1,
        'B': 3,
        'Al': 3,
        'Si': 4,
        'As': 3,
        'Se': 2,
        'Te': 2,
    }
    
    # Get the normal valence for this atom
    if atom_symbol not in normal_valence:
        return 0
    
    expected_valence = normal_valence[atom_symbol]
    
    # Special cases for common charged states
    if atom_symbol == 'N':
        if bond_count == 4:  # Quaternary ammonium
            return 1
        elif bond_count == 2:  # Anionic nitrogen
            return -1
        elif bond_count == 3:  # Normal neutral nitrogen
            return 0
    elif atom_symbol == 'O':
        if bond_count == 3:  # Oxonium
            return 1
        elif bond_count == 1:  # Alkoxide/phenoxide
            return -1
        elif bond_count == 2:  # Normal neutral oxygen
            return 0
    elif atom_symbol == 'S':
        if bond_count == 3:  # Sulfonium
            return 1
        elif bond_count == 1:  # Thiolate
            return -1
        elif bond_count in [2, 4, 6]:  # S can expand octet
            return 0
    elif atom_symbol == 'P':
        if bond_count == 4:  # Phosphonium
            return 1
        elif bond_count in [3, 5]:  # P can expand octet
            return 0
    elif atom_symbol == 'C':
        if bond_count == 3:  # Carbocation or carbanion (need more info)
            # Without lone pair info, we can't distinguish
            return 0  # Assume neutral (might be in aromatic ring)
        elif bond_count == 4:
            return 0
    elif atom_symbol in ['F', 'Cl', 'Br', 'I']:
        if bond_count == 0:  # Halide anion
            return -1
        elif bond_count == 1:  # Normal halogen
            return 0
        elif bond_count > 1:  # Hypervalent halogen (Cl, Br, I only)
            return 0
    elif atom_symbol == 'B':
        if bond_count == 4:  # Borate
            return -1
        elif bond_count == 3:
            return 0
    elif atom_symbol == 'H':
        if bond_count == 0:  # Hydride
            return -1
        elif bond_count == 1:
            return 0
    
    # Default: if bond count differs significantly from expected, calculate charge
    if bond_count != expected_valence:
        # For most organic atoms: charge = expected_valence - bond_count
        # But this is a simplification and might not always be correct
        return 0
    
    return 0


def add_charge_line_to_sdf(input_file: str, output_file: str) -> None:
    with open(input_file, 'r') as infile:
        lines: List[str] = infile.readlines()

    atom_bonds: Dict[int, int] = {}
    atom_symbols: List[str] = []
    
    # Find where the atom block starts (after the counts line)
    atom_block_start = -1
    atom_count = 0
    bond_count = 0
    
    for idx, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 10 and parts[0].isdigit() and parts[1].isdigit():
            # This is the counts line
            atom_count = int(parts[0])
            bond_count = int(parts[1])
            atom_block_start = idx + 1
            break

    
    # Parse the atom block to get atom symbols
    for idx in range(atom_block_start, atom_block_start + atom_count):
        if idx >= len(lines):
            break
        line = lines[idx]
        columns: List[str] = line.split()
        if len(columns) >= 4:
            atom_symbols.append(columns[3])
            atom_bonds[idx - atom_block_start + 1] = 0

    # Parse the bond block to count bonds for each atom
    bond_block_start = atom_block_start + atom_count
    
    for idx in range(bond_block_start, bond_block_start + bond_count):
        if idx >= len(lines):
            break
        line = lines[idx]
        if line.strip().startswith("M"):  # Stop if we hit property block
            break
        columns: List[str] = line.split()
        if len(columns) >= 3:
            try:
                atom1: int = int(columns[0])
                atom2: int = int(columns[1])
                bond_type: int = int(columns[2])
                
                # For formal charge calculation, we typically use bond order
                # Single = 1, Double = 2, Triple = 3, Aromatic = 1.5 (but stored as 4)
                bond_order = bond_type if bond_type <= 3 else 1.5
                
                if atom1 in atom_bonds:
                    atom_bonds[atom1] += bond_order
                if atom2 in atom_bonds:
                    atom_bonds[atom2] += bond_order
            except ValueError:
                continue

    # Calculate charges for each atom
    charge_positions: List[Tuple[int, int]] = []
    for atom_index, bond_count in atom_bonds.items():
        if 0 < atom_index <= len(atom_symbols):
            atom_symbol: str = atom_symbols[atom_index - 1]
            # Round bond count for comparison (in case of aromatic bonds)
            rounded_bond_count = round(bond_count)
            charge: int = calculate_charge(atom_symbol, rounded_bond_count)
            if charge != 0:
                charge_positions.append((atom_index, charge))
                print(f"Atom {atom_index} ({atom_symbol}) with {bond_count} bonds has charge {charge}")

    # Build charge line if there are charges
    charge_line = ""
    if charge_positions:
        charge_line_parts: List[str] = ["M  CHG", f"{len(charge_positions):3d}"]
        for pos, charge in charge_positions:
            charge_line_parts.append(f"{pos:4d}")
            charge_line_parts.append(f"{charge:4d}")
        charge_line = "".join(charge_line_parts) + "\n"

    # Write output file
    with open(output_file, 'w') as outfile:
        m_end_written = False
        for line in lines:
            if line.strip() == "M  END" and charge_line and not m_end_written:
                outfile.write(charge_line)
                m_end_written = True
            outfile.write(line)
            
def run_command(command, log_file_path):
    """Run a shell command, print its output, and write it to a log file."""
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
        print(result.stdout)
        
        with open(log_file_path, 'a') as log_file:
            log_file.write(result.stdout)
            
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing: {command}")
        print(e.stdout)
        print(e.stderr)
        
        with open(log_file_path, 'a') as log_file:
            log_file.write(f"An error occurred while executing: {command}\n")
            log_file.write(e.stdout)
            log_file.write(e.stderr)
        
        sys.exit(1)


def download_pdb(pdb_id, output_path):
    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, "w") as file:
                file.write(response.text)
            print(f"PDB file from RCSB saved to {output_path}")
    
    except Exception as e:
        print(e)
        url = f"https://files.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb{pdb_id.lower()}.ent.gz"
        response = requests.get(url)
        if response.status_code == 200:
            with open(f"pdb{pdb_id.lower()}.ent.gz", 'wb') as f:
                f.write(response.content)
            subprocess.run(["gunzip", f"pdb{pdb_id.lower()}.ent.gz"])
            subprocess.run(["mv", f"pdb{pdb_id.lower()}.ent", f"{output_path}"])
            print(f"[INFO] Downloaded PDB ID '{pdb_id}' to '{output_path}'")
        else:
            raise FileNotFoundError(f"Could not fetch PDB ID '{pdb_id}' from RCSB.")
    


class ProteinHOHSelect(Select):
    def __init__(self, protein_chain_id):
        self.protein_chain_id = protein_chain_id
        #self.water_chain_id = water_chain_id

    def accept_chain(self, chain):
        return chain.id == self.protein_chain_id #or chain.id == self.water_chain_id

    def accept_residue(self, residue):
        resname = residue.get_resname()
        hetfield = residue.id[0]
        chain_id = residue.get_parent().id

        is_protein = (chain_id == self.protein_chain_id and hetfield == " ")
        #is_water = (chain_id == self.water_chain_id and resname in ['HOH', 'WAT', 'SOL'])

        return is_protein


class LigandSelect(Select):
    def __init__(self, chain_id, res_id, res_name):
        self.chain_id = chain_id
        self.res_id = int(res_id)
        self.res_name = res_name

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return (
            residue.get_parent().id == self.chain_id and
            residue.get_resname().strip() == self.res_name and
            residue.id[1] == self.res_id
        )


def insert_ter_before_waters(pdb_path):
    with open(pdb_path) as f:
        pdb_lines = f.readlines()

    output_lines = []
    ter_inserted = False
    atom_serial = 1

    for i, line in enumerate(pdb_lines):
        record_type = line[:6].strip()
        is_water = (
            record_type == "HETATM" and
            line[17:20].strip() in {"HOH", "WAT", "SOL"}
        )

        if not ter_inserted and is_water:
            # Insert TER line before first water
            last_atom_line = output_lines[-1]
            resname = last_atom_line[17:20]
            chain_id = last_atom_line[21]
            resseq = last_atom_line[22:26]
            ter_line = f"TER   {atom_serial:5d}      {resname} {chain_id}{resseq}                                                  \n"
            output_lines.append(ter_line)
            atom_serial += 1
            ter_inserted = True

        if record_type in {"ATOM", "HETATM"}:
            new_line = f"{record_type:<6}{atom_serial:5d}{line[11:]}"
            output_lines.append(new_line)
            atom_serial += 1
        else:
            output_lines.append(line)

    with open(pdb_path, "w") as f:
        f.writelines(output_lines)



def export_protein_ligand_from_file(pdb_input, prot_chain_id, lig_chain_id, lig_resid, lig_resname):
    # Handle PDB ID or file
    if len(pdb_input) == 4 and not os.path.isfile(pdb_input):
        pdb_file = f"{pdb_input.upper()}.pdb"
        download_pdb(pdb_input, pdb_file)
    else:
        pdb_file = pdb_input

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file)

    base = os.path.splitext(os.path.basename(pdb_file))[0]

    # === Extract protein + water ===
    protein_file = f"prot_{prot_chain_id}.pdb"
    io_protein = PDBIO()
    io_protein.set_structure(structure)
    io_protein.save(protein_file, ProteinHOHSelect(prot_chain_id))
    insert_ter_before_waters(protein_file)

    # === Extract ligand ===
    ligand_pdb = f"lig_{lig_chain_id}_{lig_resid}_{lig_resname}.pdb"
    io_ligand = PDBIO()
    io_ligand.set_structure(structure)
    io_ligand.save(ligand_pdb, LigandSelect(lig_chain_id, lig_resid, lig_resname))

    # === Add hydrogens and write SDF ===
    h_addition = check_hydrogens_pdb(ligand_pdb)
    print(h_addition["status"])
    ligand_sdf = "prep_lig.sdf"
    if h_addition["status"] == "all hydrogens present":
        run_command(f"obabel {ligand_pdb} -O {ligand_sdf}", log_file_path="sys_prep_and_openmm_md.log")
    else:
        run_command(f"obabel {ligand_pdb} -O {ligand_sdf} -h", log_file_path="sys_prep_and_openmm_md.log")

    if len(pdb_input) == 4 and not os.path.isfile(pdb_input):
        filenames = download_sdf_files(pdb_input.upper(), lig_resname)
        for sdf_file in filenames:
            matched_sdf_file = get_matching_sdf_filename(ligand_pdb, sdf_file)
            if matched_sdf_file is not None:
                break 
    
        # Add hydrogens to the ligand
        run_command(f"obabel {matched_sdf_file} -O {ligand_sdf} -h", log_file_path="sys_prep_and_openmm_md.log")
    else:
        add_charge_line_to_sdf(ligand_sdf, ligand_sdf)

    print("\n[✔] Preparation complete. Files generated:")
    print(f"  - Protein (with water): {protein_file}")
    print(f"  - Ligand with hydrogens: {ligand_sdf}")

    return protein_file, ligand_sdf
