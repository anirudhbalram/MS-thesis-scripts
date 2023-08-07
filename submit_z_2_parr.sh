#!/bin/bash
#SBATCH --qos=narayanan-b
#SBATCH --job-name=sigma_2    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ranirudh@students.iisertirupati.ac.in     # Where to send mail
#SBATCH --nodes=1                    # Run all processes on a single node       
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --array=1-268
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=sigma_2.log   # Standard output and error log

module load conda

conda activate /blue/narayanan/a.ravishankar/py38

python projected_sigma_lists_calc_z_2_fix.py --galinfoline $SLURM_ARRAY_TASK_ID
