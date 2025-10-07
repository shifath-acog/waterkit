import requests
from Bio import PDB
from Bio.PDB import PDBParser, PPBuilder
from Bio import pairwise2
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Simulation, HBonds
from openmm import CustomExternalForce, unit, LangevinIntegrator
from typing import Set, Dict, List
from System_Prep.split_sys import export_protein_ligand_from_file, add_charge_line_to_sdf
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
from openmm import app, unit, LangevinIntegrator, Vec3
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, XTCReporter
from openmm import CustomExternalForce, System
from openmm import Platform
import subprocess
import os
from rdkit import Chem

def convert_pdb_to_sdf(input_pdb, output_sdf):
    """
    Converts a PDB file to an SDF file using Open Babel.
    """
    command = ["obabel", input_pdb, "-O", output_sdf]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
        print(f"Successfully converted {input_pdb} to {output_sdf}")
    else:
        print(f"Error converting {input_pdb} to {output_sdf}: {result.stderr.decode()}")


def download_pdb(pdb_id, output_file):
    """Downloads a PDB file from the RCSB Protein Data Bank."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_file, "w") as file:
            file.write(response.text)
        print(f"PDB file from RCSB saved to {output_file}")
    else:
        raise ValueError(f"Failed to fetch PDB file for ID {pdb_id}. HTTP status code: {response.status_code}")

def get_chain_ids_from_pdb(file_path: str) -> Set[str]:
    """Extracts all unique chain IDs from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)
    chain_ids = {chain.id for model in structure for chain in model}
    return sorted(chain_ids)

def extract_pdb_sequence(pdb_file, chain_id):
    """Extracts the sequence from the specified chain of a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ppb = PPBuilder()
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return "".join([str(pp.get_sequence()) for pp in ppb.build_peptides(chain)])
    return None

def find_missing_residues(pdb_file, chain_id, pdb_id=None):
    """Identifies missing residues by aligning the PDB sequence with the reference sequence."""
    pdb_sequence = extract_pdb_sequence(pdb_file, chain_id)
    if not pdb_sequence:
        print(f"Chain {chain_id} not found in the PDB file.")
        return None, None, None
    fasta_data = download_fasta_from_rcsb(pdb_id)
    reference_sequence = "".join(line.strip() for line in fasta_data.splitlines() if not line.startswith(">"))
    alignments = pairwise2.align.globalxs(reference_sequence, pdb_sequence, -2, -1)
    ref_aligned, pdb_aligned = alignments[0][0], alignments[0][1]
    missing_residues = {}
    ref_position = 0
    for ref_res, pdb_res in zip(ref_aligned, pdb_aligned):
        if ref_res != "-":
            ref_position += 1
        if ref_res != "-" and pdb_res == "-":
            missing_residues[ref_position] = ref_res
    return missing_residues, ref_aligned, pdb_aligned

def download_fasta_from_rcsb(pdb_id):
    """Downloads a protein FASTA sequence from the RCSB PDB database."""
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise ValueError(f"Failed to fetch sequence for PDB ID {pdb_id}. HTTP status code: {response.status_code}")

def convert_missing_residues(missing_residues: Dict[int, str], chain_index: int):
    """Converts a dictionary of missing residues into a grouped format with cumulative indices."""
    missing_residue = {}
    current_group = []
    cumulative_length = 0
    start_point = list(missing_residues.keys())[0]
    sorted_keys = sorted(missing_residues.keys())
    for i, index in enumerate(sorted_keys):
        if i == 0 or index == sorted_keys[i - 1] + 1:
            current_group.append(missing_residues[index])
        else:
            key = (chain_index, start_point - cumulative_length - 1)
            missing_residue[key] = current_group
            cumulative_length += len(current_group)
            start_point = index
            current_group = [missing_residues[index]]
    if current_group:
        key = (chain_index, start_point - cumulative_length - 1)
        missing_residue[key] = current_group
    return missing_residue

def get_atom_indices_for_missing_residues(topology, missing_residues):
    """
    Extract atom indices corresponding to the missing residues.

    Parameters:
        topology: The OpenMM topology object.
        missing_residues: Dictionary with missing residues information.

    Returns:
        List of atom indices corresponding to the missing residues.
    """
    missing_residue_indices = []
    for (chain_index, start_residue), residues in missing_residues.items():
        for residue_index, residue_name in enumerate(residues):
            # Find the chain by its index
            chain = list(topology.chains())[chain_index]
            # Calculate the residue position
            residue_position = start_residue + residue_index
            # Iterate over atoms to find the matching residue and its atoms
            for atom in topology.atoms():
                if atom.residue.chain.index == chain.index and atom.residue.index == residue_position:
                    missing_residue_indices.append(atom.index)

    return missing_residue_indices

def minimize_system(input_pdb_file: str, input_sdf_file: str, output_pdb_file: str):
    """Perform energy minimization only on missing residues using OpenMM."""
    protein_pdb = PDBFile(input_pdb_file)
    #ligand_mol = Molecule.from_file(input_sdf_file,file_format="sdf")
    rdmol = Chem.SDMolSupplier(input_sdf_file, removeHs=False)[0]
    ligand_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    ligand_mol.assign_partial_charges(partial_charge_method="am1bcc")
        # Create the Modeller instance
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)

    # Assuming ligand_mol is already defined and loaded
    lig_top = ligand_mol.to_topology()
    modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
    print(f'System has {modeller.topology.getNumAtoms()} atoms after adding ligand.')
    
    print('Preparing system')
    # Initialize a SystemGenerator using the SAGE for the ligand and tip3p for the water.
    forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'hydrogenMass': 1*unit.amu }
    system_generator = SystemGenerator(
    forcefields=["amber/ff14SB.xml","amber/tip3p_standard.xml"],
    small_molecule_forcefield="openff-2.1.0",
    molecules=[ligand_mol],
    forcefield_kwargs=forcefield_kwargs)

    # Create the system using the SystemGenerator
    system = system_generator.create_system(modeller.topology, molecules=ligand_mol)

    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picoseconds, 0.002 * unit.picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    print("Minimizing only the system...")
    simulation.minimizeEnergy()
    #minimized_positions = simulation.context.getState(getPositions=True).getPositions()
    with open(output_pdb_file, 'w') as outfile:
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)
    print(f"Minimized structure saved to {output_pdb_file}")
    #split_pdb("output_pdb_file", "chains_AB.pdb", "chain_X_UNK.pdb")

def split_pdb(prot_chain_id, input_pdb, output_pdb_1, output_pdb_2, output_sdf):
    """
    Splits a PDB file into two separate files:
    - One containing chains A and B, including HOH, WAT, SOL residues
    - Another containing chain X with resname UNK
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    io = PDB.PDBIO()
    
    # First PDB file: Chains A and B, including HOH, WAT, SOL residues
    class ChainSelectorABWithWater(PDB.Select):
        def accept_residue(self, residue):
            return residue.parent.id in [prot_chain_id] or residue.resname in ["HOH", "WAT", "SOL"]

    io.set_structure(structure)
    io.save(output_pdb_1, ChainSelectorABWithWater())

    # Second PDB file: Chain X with resname UNK
    class ChainSelectorUNK(PDB.Select):
        def accept_residue(self, residue):
            return residue.parent.id == "X" and residue.resname == "UNK"

    io.set_structure(structure)
    io.save(output_pdb_2, ChainSelectorUNK())

    print(f"Saved {output_pdb_1} with chains A and B including water residues HOH, WAT, SOL")
    print(f"Saved {output_pdb_2} with chain X (UNK)")
    # Convert PDB to SDF using Open Babel
    convert_pdb_to_sdf(output_pdb_2, output_sdf)

def minimize_missing_residues(input_pdb_file: str, input_sdf_file: str, output_pdb_file: str, missing_residue_indices: List[int]):
    """Perform energy minimization only on missing residues using OpenMM."""
    protein_pdb = PDBFile(input_pdb_file)
    #ligand_mol = Molecule.from_file(input_sdf_file,file_format="sdf")
    rdmol = Chem.SDMolSupplier(input_sdf_file, removeHs=False)[0]
    ligand_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    ligand_mol.assign_partial_charges(partial_charge_method="am1bcc")
    # Create the Modeller instance
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)

    # Assuming ligand_mol is already defined and loaded
    lig_top = ligand_mol.to_topology()

    modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
    print(f'System has {modeller.topology.getNumAtoms()} atoms after adding ligand.')


    print('Preparing system')
    # Initialize a SystemGenerator using the SAGE for the ligand and tip3p for the water.
    forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': True, 'hydrogenMass': 1*unit.amu }
    system_generator = SystemGenerator(
    forcefields=["amber/ff14SB.xml","amber/tip3p_standard.xml"],
    small_molecule_forcefield="openff-2.1.0",
    molecules=[ligand_mol],
    forcefield_kwargs=forcefield_kwargs)
   
    # Create the system using the SystemGenerator
    system = system_generator.create_system(modeller.topology, molecules=ligand_mol)
    restraint_force=  CustomExternalForce('(k/2)*periodicdistance(x, y, z, x0, y0, z0)^2')
    restraint_force.addGlobalParameter('k', 1000 * unit.kilojoules_per_mole/unit.nanometer**2)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")
    for i, pos in enumerate(protein_pdb.positions):
        if i not in missing_residue_indices:
            restraint_force.addParticle(i, pos.value_in_unit(unit.nanometer))
    system.addForce(restraint_force)
    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picoseconds, 0.002 * unit.picoseconds)
    # Check and print the platform in use
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'single', 'CudaDeviceIndex': '0'}
    simulation = Simulation(modeller.topology, system, integrator, platform=platform)
    context = simulation.context
    context.setPositions(modeller.positions)

    # Check and print the platform in use
    platform_name = simulation.context.getPlatform().getName()
    print(f"Running on platform: {platform_name}")


    print('Minimising ...')
    simulation.minimizeEnergy()

    
    with open(output_pdb_file, 'w') as outfile:
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file=outfile, keepIds=True)
    print(f"Minimized structure saved to {output_pdb_file}")


def system_prep(pdb_input, prot_chain_id, lig_chain_id, lig_resid, lig_resname):
    pdb_file, ligand_sdf_file = export_protein_ligand_from_file(pdb_input, prot_chain_id, lig_chain_id, lig_resid, lig_resname)
    fixed_pdb_file = "prep_comp.pdb"
    minimized_pdb_file = "prep_min_comp.pdb"
    chain_x_file="prep_min_lig.pdb"
    if len(pdb_input) == 4 and not os.path.isfile(pdb_input):
        pdb_id = pdb_input
    
    else:
        pdb_id = None

    chain_ids = get_chain_ids_from_pdb(pdb_file)
    missing_residues_pdbfixer_all = {}
    chain_index = 0 
    for chain_id in chain_ids:
        pdb_sequence = extract_pdb_sequence(pdb_file, chain_id)
        missing, ref_aligned, pdb_aligned = find_missing_residues(pdb_file, chain_id, pdb_id)
        amino_acid_map = {
            "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
            "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
            "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
            "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"
        }
        
        print(f"For {pdb_id} for chain {chain_id}")
        print("-----Reference")
        print(ref_aligned)
        print("-----pdb_sequence")
        print(pdb_aligned)
        
        seq_len = len(pdb_sequence) if pdb_sequence else 0
        
        if missing:
            missing_residues = {pos: amino_acid_map[res] for pos, res in missing.items()}
            missing_residues_pdbfixer = convert_missing_residues(missing_residues, chain_index)
            
            missing_residues_pdbfixer.pop((chain_index, 0), None)
            missing_residues_pdbfixer.pop((chain_index, seq_len), None)
            missing_residues_pdbfixer_all.update(missing_residues_pdbfixer)
        
        chain_index += 1

    fixer = PDBFixer(filename=pdb_file)
    #numChains = len(list(fixer.topology.chains()))
    #fixer.removeChains(range(1, numChains))
    fixer.findMissingResidues()

    if missing_residues_pdbfixer_all:
        fixer.missingResidues = missing_residues_pdbfixer_all
        print(f"Missing Residues Added: {fixer.missingResidues}")

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    fixer.removeHeterogens(keepWater=True)


    chain_ab_file="prep_min_prot_wat.pdb"
    chain_x_sdf="prep_min_lig.sdf"
    with open(fixed_pdb_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f,keepIds=True)

    if missing_residues_pdbfixer_all:
        missing_residue_indices = get_atom_indices_for_missing_residues(fixer.topology, fixer.missingResidues)
        if missing_residue_indices:
            minimize_missing_residues(fixed_pdb_file, ligand_sdf_file, minimized_pdb_file, missing_residue_indices)
            split_pdb(prot_chain_id, minimized_pdb_file, chain_ab_file, chain_x_file, chain_x_sdf)
            print("Missing residue indices found. Minimization of missing residues")
            minimize_system(chain_ab_file, chain_x_sdf, minimized_pdb_file)
            split_pdb(prot_chain_id, minimized_pdb_file, chain_ab_file, chain_x_file,chain_x_sdf)
            add_charge_line_to_sdf(chain_x_sdf, chain_x_sdf)
            print("Only minimisation of the system.")
        else:
            print("No missing residue indices found. Skipping minimization.")
    else:
        minimize_system(fixed_pdb_file, ligand_sdf_file, minimized_pdb_file)
        split_pdb(prot_chain_id, minimized_pdb_file, chain_ab_file, chain_x_file,chain_x_sdf)
        add_charge_line_to_sdf(chain_x_sdf, chain_x_sdf)
        print("No missing residues detected. Only minimisation of the system.")


