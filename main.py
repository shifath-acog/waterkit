import yaml
from System_Prep.split_sys import run_command
import os
from wk_utils import process_pdb_and_run_tleap
import shutil


def load_yaml(yml_file):
    with open(yml_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

config = load_yaml('setting.yml')
pdb_file = config.get("pdb_file",None)
pdb_id = config.get("pdb_id",None)
prot_chain_id = config.get("prot_chain_id",None)
lig_chain_id = config.get("lig_chain_id",None)
lig_resid = config.get("lig_resid",None)
lig_resname = config.get("lig_resname",None)
lig_charge = config.get("lig_charge",0)
wat_sites = config.get("wat_sites",None)
wat_chain_id = config.get("wat_chain_id","A")
wat_sel_dist = config.get("wat_sel_dist",5.0)
wat_model =  config.get("wat_model","tip3p")
wk_frames = config.get("wk_frames",100)
padding = config.get("padding",6.0)
n_cores = config.get("n_cores",16)
gist_res = config.get("gist_res",0.5)
output_folder = config.get("output_folder",None)
if output_folder:
    os.makedirs(output_folder, exist_ok=True)
    os.chdir(output_folder)
    shutil.copy("../gistpp.sh","./")
else:
    raise ValueError("Output folder ('md') is not set in the YAML file.")

if pdb_id:
    input_pdb = pdb_id
else:
    input_pdb = pdb_file

process_pdb_and_run_tleap(input_pdb, prot_chain_id=prot_chain_id, lig_chain_id=lig_chain_id, lig_resid=lig_resid, 
                            lig_resname=lig_resname, lig_charge=lig_charge, padding=padding, water_model=wat_model, 
                            wk_frames=wk_frames, n_cores=n_cores, gist_res=gist_res, skip_h_add=False,
                            wat_chain_id=wat_chain_id, wat_sel_dist=wat_sel_dist, wat_sites=wat_sites)

shutil.copy("../setting.yml","./")