#!/bin/bash
#SBATCH --qos=narayanan-b
#SBATCH --job-name=half_radii_2    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ranirudh@students.iisertirupati.ac.in     # Where to send mail
#SBATCH --nodes=1                    # Run all processes on a single node       
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --array=1-268
#SBATCH --mem=80gb                     # Job memory request
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=half_radii_2.log   # Standard output and error log

module load conda

conda activate /blue/narayanan/a.ravishankar/py38

python half_radii_z_2_fix.py --galinfoline $SLURM_ARRAY_TASK_ID
