#!/bin/bash
#SBATCH -N 1
#SBATCH -o run.out
#SBATCH -e run.err
#SBATCH -w atmnode004
#SBATCH --mem-per-cpu=2000M
#SBATCH --partition=priority-rp

python run.py
