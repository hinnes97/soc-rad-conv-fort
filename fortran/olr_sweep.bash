#!/bin/bash
#SBATCH --output=final_mstar.out
#SBATCH --error=final_mstar.err
#SBATCH --partition=priority-rp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -w atmnode017
#SBATCH --mem-per-cpu=1000

python -u olr_sweep.py

