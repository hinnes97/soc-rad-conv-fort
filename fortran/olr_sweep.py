import numpy as np
import os
import xarray as xr
S_vec = np.arange(50,200,10.)


directory = 'inhibited_output'
os.system('mkdir '+directory)

for i,S in enumerate(S_vec):
    fname = f'olr_{i}'
    input_file = f'olr_{i-1}.nc'
    command = 'python run.py ' + directory + ' ' + fname + ' ' + str(S) + ' ' + input_file
    
    os.system(command)
    test = xr.open_dataset(directory+'/'+fname+'.nc')
    if (np.absolute((S - test.olr)/S) > 0.1):
        test.close()
        os.system('rm '+directory+'/'+fname+'.nc')
        i_end = i-1
        break
    i_end = i

np.savetxt(directory+'/S_vec.txt', S_vec[:i_end+1])

os.system('python merge_data.py '+ directory + ' 60 61')
