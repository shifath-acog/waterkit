import uuid
import json
import fcntl
import yaml
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Optional, List

from fastapi import FastAPI, HTTPException, File, UploadFile
from fastapi.responses import FileResponse
from pydantic import BaseModel

# --- Configuration ---
STATUS_FILE = Path("/home/waterkit/WATERKIT/task-status.json")
OUTPUTS_DIR = Path("/home/waterkit/WATERKIT/outputs")
UPLOADS_DIR = Path("/home/waterkit/WATERKIT/uploads")
PROJECT_ROOT = Path("/home/waterkit/WATERKIT")


# --- Pydantic Models ---
# Defines the structure of the API request body
class PipelineSettings(BaseModel):
    pdb_file: Optional[str] = None
    pdb_id: Optional[str] = None
    prot_chain_id: str
    lig_chain_id: str
    lig_resid: int
    lig_resname: str
    lig_charge: int = 0
    wat_chain_id: str = "A"
    wat_sel_dist: float = 5.0
    wat_sites: Optional[List[str]] = None
    wat_model: str = "tip3p"
    wk_frames: int = 100
    padding: float = 6.0
    n_cores: int = 16
    gist_res: float = 0.5

    class Config:
        extra = 'ignore' # Ignore extra fields from request

# --- FastAPI App ---
app = FastAPI()

@app.get("/")
def read_root():
    return {"message": "WaterKit API is running"}

@app.post("/api/upload")
def upload_file(file: UploadFile = File(...)):
    """
    Accepts a file upload and saves it to the 'uploads' directory.
    Returns the path to the saved file.
    """
    try:
        # Ensure the uploads directory exists
        UPLOADS_DIR.mkdir(exist_ok=True)

        # Sanitize filename to prevent directory traversal attacks
        filename = Path(file.filename).name
        if not filename:
            raise HTTPException(status_code=400, detail="Invalid filename.")
            
        dest_path = UPLOADS_DIR / filename

        # Save the file
        with dest_path.open("wb") as buffer:
            shutil.copyfileobj(file.file, buffer)

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to upload file: {e}")
    finally:
        file.file.close()

    return {"file_path": str(dest_path)}


@app.post("/api/tasks")
def create_task(settings: PipelineSettings):
    """
    Creates and starts a new pipeline task.
    """
    if not settings.pdb_file and not settings.pdb_id:
        raise HTTPException(status_code=400, detail="Either 'pdb_file' or 'pdb_id' must be provided.")

    task_id = str(uuid.uuid4())
    task_dir = OUTPUTS_DIR / task_id
    
    # 1. Create a dedicated directory for the task
    try:
        task_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise HTTPException(status_code=500, detail=f"Failed to create task directory: {e}")

    # 2. Generate the setting.yml file for this specific task
    settings_dict = settings.dict()
    # The output_folder for the script is the task directory itself.
    settings_dict['output_folder'] = str(task_dir)

    # Ensure input file paths are absolute within the container
    if settings.pdb_file and not settings.pdb_file.startswith('/'):
        settings_dict['pdb_file'] = f"/home/waterkit/WATERKIT/{settings.pdb_file}"

    setting_yml_path = task_dir / "setting.yml"
    with open(setting_yml_path, 'w') as f:
        yaml.dump(settings_dict, f)

    # Copy gistpp.sh to the task directory
    shutil.copy(PROJECT_ROOT / "gistpp.sh", task_dir)

    # 3. Log the task in the central status file (with file locking)
    new_task_info = {
        "task_id": task_id,
        "status": "QUEUED",
        "submitted_at": datetime.utcnow().isoformat(),
        "task_dir": str(task_dir),
        "params": settings.dict()
    }

    try:
        # Open the file in read-write mode, creating it if it doesn't exist
        with open(STATUS_FILE, "a+") as f:
            f.seek(0) # Move to the beginning to read
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                data = {}
            
            data[task_id] = new_task_info
            
            f.seek(0) # Move back to the beginning to write
            f.truncate() # Clear the file before writing
            json.dump(data, f, indent=4)
            
            fcntl.flock(f, fcntl.LOCK_UN)

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to write to status file: {e}")

    # 4. Launch the main.py script as a background process
    command = ["/miniconda/envs/waterkit/bin/python", "/home/waterkit/WATERKIT/main.py", "--task-id", task_id]
    # Get current environment variables
    import os
    env = os.environ.copy()
    # Prepend the conda environment's bin directory to PATH
    env["PATH"] = "/miniconda/envs/waterkit/bin:" + env["PATH"]

    subprocess.Popen(command, env=env)

    return {"task_id": task_id, "message": "Task successfully queued."}

@app.get("/api/tasks/{task_id}/status")
def get_task_status(task_id: str):
    """
    Retrieves the current status of a specific task.
    """
    try:
        with open(STATUS_FILE, "r") as f:
            fcntl.flock(f, fcntl.LOCK_SH) # Use shared lock for reading
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                data = {} # If file is empty or malformed, treat as empty dict
            fcntl.flock(f, fcntl.LOCK_UN)

        task_info = data.get(task_id)
        if not task_info:
            raise HTTPException(status_code=404, detail=f"Task with ID {task_id} not found.")
        
        return task_info

    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Task status file not found.")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to read task status: {e}")

@app.get("/api/tasks/{task_id}/results")
def get_task_results(task_id: str):
    """
    Retrieves a list of output files for a completed task.
    """
    try:
        with open(STATUS_FILE, "r") as f:
            fcntl.flock(f, fcntl.LOCK_SH) # Use shared lock for reading
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                data = {}
            fcntl.flock(f, fcntl.LOCK_UN)

        task_info = data.get(task_id)
        if not task_info:
            raise HTTPException(status_code=404, detail=f"Task with ID {task_id} not found.")
        
        if task_info["status"] != "COMPLETED":
            raise HTTPException(status_code=409, detail=f"Task {task_id} is not completed yet. Current status: {task_info['status']}.")

        task_dir = Path(task_info["task_dir"])
        if not task_dir.is_dir():
            raise HTTPException(status_code=404, detail=f"Output directory for task {task_id} not found.")

        # List files in the task directory, excluding specific internal files
        excluded_files = ["setting.yml", "gistpp.sh", "logs"]
        result_files = []
        for item in task_dir.iterdir():
            if item.name not in excluded_files and item.is_file():
                result_files.append(item.name)
            elif item.name not in excluded_files and item.is_dir():
                # Optionally, list contents of subdirectories or just include directory name
                result_files.append(item.name + "/") # Indicate it's a directory

        return {"task_id": task_id, "status": task_info["status"], "results": sorted(result_files)}

    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Task status file not found.")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to retrieve task results: {e}")

@app.get("/api/tasks/{task_id}/files/{filename}")
def get_task_file(task_id: str, filename: str):
    """
    Retrieves the content of a specific output file for a completed task.
    """
    try:
        with open(STATUS_FILE, "r") as f:
            fcntl.flock(f, fcntl.LOCK_SH) # Use shared lock for reading
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                data = {}
            fcntl.flock(f, fcntl.LOCK_UN)

        task_info = data.get(task_id)
        if not task_info:
            raise HTTPException(status_code=404, detail=f"Task with ID {task_id} not found.")

        if task_info["status"] != "COMPLETED":
            raise HTTPException(status_code=409, detail=f"Task {task_id} is not completed yet. Current status: {task_info['status']}.")

        task_dir = Path(task_info["task_dir"])
        file_path = task_dir / filename

        # Ensure the file exists and is within the task directory
        if not file_path.is_file() or not file_path.is_relative_to(task_dir):
            raise HTTPException(status_code=404, detail=f"File {filename} not found for task {task_id}.")

        return FileResponse(file_path, media_type="application/octet-stream", filename=filename)

    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Task status file not found.")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to retrieve file content: {e}")