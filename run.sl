#!/bin/bash -e
#SBATCH --job-name=phyloParrot    # job name (shows up in queue)
#SBATCH --time=00-01:00:00        # Walltime (DD-HH:MM:SS)
#SBATCH --mem=20000               # total memory in MB
#SBATCH --cpus-per-task=4         # 4 CPUs
#SBATCH --account=uoa03415        # Project code

# load R
module load R

# run script
Rscript run.R
