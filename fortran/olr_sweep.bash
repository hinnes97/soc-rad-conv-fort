#!/bin/bash
#SBATCH --output=mstar_tweak.out
#SBATCH --error=mstar_tweak.err
#SBATCH --partition=priority-rp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -w atmnode001
#SBATCH --mem-per-cpu=1000

python -u olr_sweep.py

