#!/bin/bash
#SBATCH --output=final_mstar_10bar.out
#SBATCH --error=final_mstar_10bar.err
#SBATCH --partition=priority-rp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -w atmnode004
#SBATCH --mem-per-cpu=1000

stdbuf -o0 -e0 python -u olr_sweep.py

