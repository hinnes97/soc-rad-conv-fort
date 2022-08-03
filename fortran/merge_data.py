import xarray as xr
import numpy as np
import argparse
import os
import re
import glob

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

    
parser = argparse.ArgumentParser(description='input arguments')
parser.add_argument('direc', type=str, help='file directory')
parser.add_argument('nf', type=int, help='Number of full levels')
parser.add_argument('ne', type=int, help = 'Number of edge levels')
args =parser.parse_args()

# Load instellations
S_vec = np.loadtxt(args.direc+'/S_vec.txt')
#S_vec = np.linspace(10,20,2)
files = sorted(glob.glob(args.direc+'/*.nc'), key=numericalSort)
dss = []
for infile in files:
    dss.append(xr.open_dataset(infile))

output = xr.concat(dss, xr.DataArray(S_vec, coords={'S':S_vec}))
os.system('mkdir ' + args.direc + '/merged')    
output.to_netcdf(args.direc+'/merged/merged_data.nc')
    

