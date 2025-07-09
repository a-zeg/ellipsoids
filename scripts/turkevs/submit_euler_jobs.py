import os
import sys
import subprocess
import argparse

sys.path.append(os.path.abspath('.'))

from ellipsoids.turkevs.config import TURKEVS_DATA_DIR



def create_sbatch_script(folder_path, job_name):
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
python ~/ellipsoids/scripts/turkevs/calculate_ellipsoids_barcodes_turkevs.py --folder "{folder_path}"
"""



def submit_job(folder_path, dry_run=False):
    job_name = f"{os.path.basename(folder_path)}"
    script_content = create_sbatch_script(folder_path, job_name)

    if dry_run:
        print(f"\n[SIMULATION] Would submit job for: {folder_path}")
        print("=" * 60)
        print(script_content)
        print("=" * 60)
        return

    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".sh") as f:
        f.write(script_content)
        script_path = f.name

    try:
        subprocess.run(["sbatch", script_path], check=True)
        print(f"Submitted: {folder_path}")

    finally:
        os.remove(script_path)



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--dry_run", action="store_true", help="Print sbatch script instead of submitting it")
    args = parser.parse_args()


    if not os.path.isdir(TURKEVS_DATA_DIR):
        raise FileNotFoundError(f"Data directory does not exist: {TURKEVS_DATA_DIR}")

    for name in os.listdir(TURKEVS_DATA_DIR):
        path = os.path.join(TURKEVS_DATA_DIR, name)
        if os.path.isdir(path):
            submit_job(path, dry_run=args.dry_run)



if __name__ == "__main__":
    main()
