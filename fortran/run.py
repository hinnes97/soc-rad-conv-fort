import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('direc', type=str)
parser.add_argument('f_name', type=str)
parser.add_argument('S', type = float)
parser.add_argument('input_file',type=str)
args = parser.parse_args()

# Experiment name
#exp_name = "runaway_test"
# Set output directory
#output_dir = 'output/runaway'
output_dir = args.direc

# Ensure output directory exists
if (os.path.isdir(output_dir)==False):
    os.system('mkdir -p '+ output_dir)




# Build executable
err = os.system('./build -r SOC')

if (err != 0):
    exit('ERROR: Fortran compilation failed')


# Force me to write readme for each experiment
exp_name = args.f_name #input('Enter experiment name:\n')
readme = ' '#input('Write description of experiment:\n')

output_file = exp_name+'.nc'
# Set up log file
log_file = exp_name+'.log'

readme_file = exp_name+'.readme'

with open(readme_file, 'w') as f:
    f.write(readme+'\n')

# Use sed to change the namelist for the io options
# For physical parameters, edit the namelist directly
sed_1 = 'sed -i "s|^\s*output_file.*|  output_file = \'' + output_dir + '/' + output_file + '\'|g" input.nml'
sed_2 = 'sed -i "s|^\s*Finc.*|  Finc = ' + str(args.S) + '|g" input.nml'
sed_3 = 'sed -i "s|^\s*stellar_constant.*|  stellar_constant = ' + str(args.S) + '|g" input_soc.nml'
os.system(sed_1)
os.system(sed_2)
os.system(sed_3)

if os.path.isfile(args.direc+'/'+args.input_file):
    sed_3 = 'sed -i "s|^\s*init_from_file.*|  init_from_file = .true.|g" input.nml'
    sed_4 = 'sed -i "s|^\s*input_file.*| input_file =\''+args.direc+'/' + args.input_file+'\'|g" input.nml'
    
else:
    sed_3 = 'sed -i "s|^\s*init_from_file.*|  init_from_file = .false.|g" input.nml'
    sed_4 = ''
os.system(sed_3)
os.system(sed_4)

os.system('cat input.nml')
# Save namelist to log file
os.system('cat input.nml > ' + log_file)

# Run program and save log to logfile, whilst also printing to stdout
#err = os.system('./main.exe 2>&1 | tee -a ' + log_file)

err = os.system('./main.exe')

if (err != 0):
    exit('ERROR: Running code failed')

# Move log and readme
os.system('mv '+ log_file + ' ' + output_dir)
os.system('mv ' + readme_file + ' ' + output_dir)
