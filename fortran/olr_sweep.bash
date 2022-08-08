#!/bin/bash
#SBATCH --output=inhibit_runaway.out
#SBATCH --error=inhibit_runaway.err
#SBATCH --partition=priority-rp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=1000

python olr_sweep.py
