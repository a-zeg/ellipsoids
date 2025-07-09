import os
import sys
import subprocess
import argparse
import tempfile
import math

sys.path.append(os.path.abspath('.'))

from ellipsoids.turkevs.turkevs_utils import read_turkevs_datasets
from ellipsoids.turkevs.config import REL_DATASETS_DIR, TURKEVS_DATA_DIR



def create_sbatch_script(folder_path, job_name, datasets_start, datasets_end):
    return f"""#!/bin/bash
#SBATCH -n 4
#SBATCH --time=28:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH --tmp=4000
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err
#SBATCH --mail-type=END,FAIL

hostname
python ~/ellipsoids/scripts/turkevs/1-calculate_ellipsoid_barcodes.py --folder "{folder_path}" --datasets_start {datasets_start} --datasets_end {datasets_end}
"""



def submit_job(folder_path, datasets_start, datasets_end, dry_run=False):
    job_name = f"{os.path.basename(folder_path)}_datasets={datasets_start}-{datasets_end}"
    script_content = create_sbatch_script(folder_path, job_name, datasets_start, datasets_end)

    if dry_run:
        print(f"\n[SIMULATION] Would submit job for: {folder_path}")
        print("=" * 60)
        print(script_content)
        print("=" * 60)
        return

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".sh") as f:
        f.write(script_content)
        script_path = f.name

    try:
        subprocess.run(["sbatch", script_path], check=True)
        print(f"Submitted: {folder_path} datasets {datasets_start} to {datasets_end}")

    finally:
        os.remove(script_path)



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--dry_run", action="store_true", help="Print sbatch script instead of submitting it")
    args = parser.parse_args()

    # folder = "turkevs_npc=100_np=300_s=0"
    folder = "test"
    n_datasets_chunks = 10

    if not os.path.isdir(TURKEVS_DATA_DIR):
        raise FileNotFoundError(f"Data directory does not exist: {TURKEVS_DATA_DIR}")
    path = os.path.join(TURKEVS_DATA_DIR, folder)

    if os.path.isdir(path):
        datasets_folder = os.path.join(path, REL_DATASETS_DIR)
        datasets_file = [os.path.join(datasets_folder, f) for f in os.listdir(datasets_folder) if "datasets" in f][0]
        datasets = read_turkevs_datasets(datasets_file)
        n_datasets = len(datasets)

        datasets_chunk_size = math.ceil(n_datasets / n_datasets_chunks)

        for i in range(n_datasets_chunks):
            start = i*datasets_chunk_size
            end = min(start + datasets_chunk_size, n_datasets)
            if start >= end:
                continue
            submit_job(path, start, end, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
