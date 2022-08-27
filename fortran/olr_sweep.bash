#!/bin/bash
#SBATCH --output=final_mstar_2.out
#SBATCH --error=final_mstar_2.err
#SBATCH --partition=priority-rp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -w atmnode022
#SBATCH --mem-per-cpu=1000

python -u olr_sweep.py

