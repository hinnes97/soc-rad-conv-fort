import numpy as np
import os
import argparse
import subprocess
import sys
import f90nml

# Build executable, change -r SOC for other radiation scheme as required

# Only want to invoke cmake if there is a new file in the src 
os.system('find src -type f \( -path "*/socrates/src*" -prune -o -print \) > ./build/srclist_new.txt')
err = os.system('cmp -s ./build/srclist_new.txt ./build/srclist.txt')
if (err != 0):
    # Source list has changed, re-invoke cmake
    err = os.system('cd build && cmake ..')
    if (err != 0):
        exit('ERROR: CMake failed')

os.system('mv ./build/srclist_new.txt ./build/srclist.txt')
err = os.system('cd build && make')
if (err!=0):
    exit('Error: Fortran compilation failed')
    

nml = f90nml.read('input.nml')
exp_name = nml['io_nml']['output_file'].split('/')[-1].split('.')[0]
output_dir = '/'.join(nml['io_nml']['output_file'].split('/')[:-1])

os.system('mkdir -p ' + output_dir)
readme = ' '#input('Write description of experiment:\n')

# Set up log file
log_file = exp_name+'.log'
readme_file = exp_name+'.readme'


readme = 'Some message for this run'
with open(readme_file, 'w') as f:
    f.write(readme+'\n')

os.system('cat input.nml')
# Save namelist to log file
os.system('cat input.nml > ' + log_file)

# Run program
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

# Plot scripts
infile = '../'+nml['io_nml']['output_file']
outfile = '../'+'/'.join(nml['io_nml']['output_file'].split('/')[:-1]) + '/ptqf.pdf'
os.system(f"cd plot_scripts && python plot_flux.py {infile} {outfile}")
sys.exit(0)

