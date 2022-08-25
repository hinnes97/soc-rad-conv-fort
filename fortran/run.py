import numpy as np
import os
import argparse
import subprocess
import sys
import f90nml

# Build executable
err = os.system('./build -r SOC')

if (err != 0):
    exit('ERROR: Fortran compilation failed')

nml = f90nml.read('input.nml')
exp_name = nml['io_nml']['output_file'].split('/')[-1].split('.')[0]
output_dir = nml['io_nml']['output_file'].split('/')[0]
readme = ' '#input('Write description of experiment:\n')

# Set up log file
log_file = exp_name+'.log'
readme_file = exp_name+'.readme'

with open(readme_file, 'w') as f:
    f.write(readme+'\n')

os.system('cat input.nml')
# Save namelist to log file
os.system('cat input.nml > ' + log_file)

# Run program and save log to logfile, whilst also printing to stdout
#err = os.system('./main.exe 2>&1 | tee -a ' + log_file)

with open(log_file, 'a') as f:
    proc = subprocess.Popen(['./main.exe'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
    for line in proc.stdout:
        sys.stdout.write(line)
        f.write(line)
    proc.wait()

if proc.returncode == 99:
    # Exit with equivalent code for olr_sweep.py
    sys.exit(99)
elif proc.returncode==100:
    sys.exit(100)
elif proc.returncode==101:
    sys.exit(101)

# Move log and readme
os.system('mv '+ log_file + ' ' + output_dir)
os.system('mv ' + readme_file + ' ' + output_dir)
sys.exit(0)
