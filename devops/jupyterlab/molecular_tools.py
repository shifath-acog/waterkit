#!/usr/bin/env python3

import os
import glob
import subprocess
import argparse

def process_pdb(input_pdb, ligand_name):
    print(f"\033[1m\033[34mProcessing PDB file: \033[92m{input_pdb}\033[0m")
    with open(input_pdb, 'r') as file:
        lines = file.readlines()
    
    print(f"\033[1m\033[34mFiltering ligand lines for ligand: \033[92m{ligand_name}\033[0m")
    filtered_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM', 'TER'))]
    ligand_lines = [line for line in lines if line.startswith('HETATM') and ligand_name in line]
    
    temp_ligand_pdb = 'temp_ligand.pdb'
    with open(temp_ligand_pdb, 'w') as file:
        file.writelines(ligand_lines)
    
    print("\033[1m\033[34mAdding hydrogens to the ligand\033[0m")
    temp_with_hydrogens = 'temp_with_hydrogens.pdb'
    subprocess.run(f"obabel {temp_ligand_pdb} -O {temp_with_hydrogens} -h",
                   shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.remove(temp_ligand_pdb)
    
    with open(temp_with_hydrogens, 'r') as file:
        lines = file.readlines()
    os.remove(temp_with_hydrogens)
    
    filtered_ligand_lines = [line for line in lines if line.startswith('HETATM') and ligand_name in line]
    ligand_with_hydrogens = 'Ligand_wH.pdb'
    with open(ligand_with_hydrogens, 'w') as file:
        file.writelines(filtered_ligand_lines)
    
    start_index = None
    end_index = None
    for i, line in enumerate(filtered_lines):
        if line.startswith('HETATM') and ligand_name in line:
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
    
    print(f"\033[1m\033[34mFinal PDB with hydrogens saved as: \033[92m{output_pdb}\033[0m")
    print(f"\033[1m\033[34mLigand PDB with hydrogens saved as: \033[92m{ligand_with_hydrogens}\033[0m")
    return output_pdb, ligand_with_hydrogens

def run_tleap(input_pdb):
    print(f"\033[1m\033[34mStarting TLEaP processing for: \033[92m{input_pdb}\033[0m")
    ligand_name = os.path.splitext(input_pdb)[0]
    
    print("\033[1m\033[34mRunning antechamber to generate MOL2 file\033[0m")
    subprocess.run(f"antechamber -i {input_pdb} -fi pdb -o {ligand_name}.mol2 -fo mol2 -c bcc -nc 0 -m 1 -at gaff2",
                   shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    print("\033[1m\033[34mRunning parmchk2 to generate FRCMOD file\033[0m")
    subprocess.run(f"parmchk2 -i {ligand_name}.mol2 -f mol2 -o {ligand_name}.frcmod",
                   shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    print("\033[1m\033[34mRunning TLEaP to prepare AMBER input files\033[0m")
    tleap_commands = (f"source leaprc.gaff; "
                      f"loadamberparams {ligand_name}.frcmod; "
                      f"BFS = loadmol2 {ligand_name}.mol2; "
                      f"check BFS; "
                      f"setBox BFS vdw 6.0; "
                      f"saveamberparm BFS {ligand_name}.prmtop {ligand_name}.inpcrd; "
                      f"saveoff BFS {ligand_name}.lib; "
                      f"savepdb BFS {ligand_name}_box.pdb; "
                      f"quit")
    
    subprocess.run(f"echo '{tleap_commands}' | tleap > leap.log",
                   shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    print("\033[1m\033[34mCleaning up temporary and intermediate files\033[0m")
    for pattern in ["ANTECHAMBER*", "sqm*", "*.inpcrd", "*.prmtop", "*.log", "*.INF"]:
        for file_path in glob.glob(pattern):
            os.remove(file_path)
    
    print(f"\033[1m\033[34mFinal PDB saved as: \033[92m{ligand_name}_box.pdb\033[0m")
    print(f"\033[1m\033[34mFinal LIB file saved as: \033[92m{ligand_name}.lib\033[0m")
    return f"{ligand_name}_box.pdb", f"{ligand_name}.lib"

def main():
    parser = argparse.ArgumentParser(description="Molecular tools for processing PDB files and running TLEaP.")
    parser.add_argument("function", type=str, choices=["process_pdb", "run_tleap"], help="Function to execute")
    parser.add_argument("input_pdb", type=str, help="Input PDB file")
    parser.add_argument("--ligand_name", type=str, help="Ligand name (required for process_pdb)", default=None)
    args = parser.parse_args()
    
    if args.function == "process_pdb":
        if not args.ligand_name:
            print("\033[1m\033[91mError: --ligand_name is required for process_pdb\033[0m")
            return
        process_pdb(args.input_pdb, args.ligand_name)
    elif args.function == "run_tleap":
        run_tleap(args.input_pdb)

if __name__ == "__main__":
    main()
