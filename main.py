import yaml
from System_Prep.split_sys import run_command
import os
from wk_utils import process_pdb_and_run_tleap
import shutil
import argparse
import json
import fcntl
import traceback
from datetime import datetime
from pathlib import Path

STATUS_FILE = "/home/waterkit/WATERKIT/task-status.json"

def update_status(task_id: str, status: str, error_message: str = None):
    """Safely updates the status of a task in the central status file."""
    try:
        with open(STATUS_FILE, "r+") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                # This case should ideally not happen if API creates the file first
                data = {}

            if task_id in data:
                task_data = data[task_id]
                task_data["status"] = status
                task_data["updated_at"] = datetime.utcnow().isoformat()
                if status == "RUNNING":
                    task_data["started_at"] = datetime.utcnow().isoformat()
                elif status in ["COMPLETED", "FAILED"]:
                    task_data["ended_at"] = datetime.utcnow().isoformat()
                
                if error_message:
                    task_data["error"] = error_message

                f.seek(0)
                f.truncate()
                json.dump(data, f, indent=4)
    except FileNotFoundError:
        # Handle case where status file doesn't exist, though API should create it
        pass  # Or log an error
    except Exception as e:
        # Log this exception, as it's an error in the status update itself
        print(f"Critical error updating status file: {e}")
    finally:
        # Ensure lock is released if file was opened
        try:
            fcntl.flock(f, fcntl.LOCK_UN)
        except Exception:
            pass

def load_yaml(yml_file):
    with open(yml_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def main_pipeline(task_id: str):
    """The core scientific pipeline logic."""
    print(f"Current working directory: {os.getcwd()}")

    # Construct the task-specific output folder path
    OUTPUTS_DIR = Path("/home/waterkit/WATERKIT/outputs")
    output_folder = OUTPUTS_DIR / task_id

    if output_folder.exists() and output_folder.is_dir():
        os.chdir(output_folder)
        print(f"Changed working directory to: {os.getcwd()}")  # Debug print
    else:
        raise ValueError(f"Task output folder not found: {output_folder}")

    # Now load the task-specific setting.yml
    config = load_yaml('setting.yml')
    pdb_file = config.get("pdb_file", None)
    pdb_id = config.get("pdb_id", None)
    prot_chain_id = config.get("prot_chain_id", None)
    lig_chain_id = config.get("lig_chain_id", None)
    lig_resid = config.get("lig_resid", None)
    lig_resname = config.get("lig_resname", None)
    lig_charge = config.get("lig_charge", 0)
    wat_sites = config.get("wat_sites", None)
    wat_chain_id = config.get("wat_chain_id", "A")
    wat_sel_dist = config.get("wat_sel_dist", 5.0)
    wat_model = config.get("wat_model", "tip3p")
    wk_frames = config.get("wk_frames", 100)
    padding = config.get("padding", 6.0)
    n_cores = config.get("n_cores", 16)
    gist_res = config.get("gist_res", 0.5)
    if pdb_id:
        input_pdb = pdb_id
    else:
        input_pdb = pdb_file

    process_pdb_and_run_tleap(
        input_pdb,
        prot_chain_id=prot_chain_id,
        lig_chain_id=lig_chain_id,
        lig_resid=lig_resid,
        lig_resname=lig_resname,
        lig_charge=lig_charge,
        padding=padding,
        water_model=wat_model,
        wk_frames=wk_frames,
        n_cores=n_cores,
        gist_res=gist_res,
        skip_h_add=False,
        wat_chain_id=wat_chain_id,
        wat_sel_dist=wat_sel_dist,
        wat_sites=wat_sites
    )
    # The final copy is removed as setting.yml is already in the correct place.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the WaterKit pipeline for a specific task.")
    parser.add_argument("--task-id", required=True, help="The unique identifier for this task.")
    args = parser.parse_args()

    task_id = args.task_id

    try:
        update_status(task_id, "RUNNING")
        main_pipeline(task_id)
        update_status(task_id, "COMPLETED")
    except Exception as e:
        print(f"Pipeline failed for task {task_id}: {e}")
        error_details = traceback.format_exc()
        update_status(task_id, "FAILED", error_message=error_details)