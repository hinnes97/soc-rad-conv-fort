import numpy as np
import os
import xarray as xr
import subprocess
import f90nml
import sys

S_vec = np.arange(50.,150.,5.)
S_vec = np.r_[np.arange(50.,130.,5.)]

directory = 'final_mstar_test'
spectral_file_lw = 'spectral_files/sp_lw_30_sub_nep_final'
spectral_file_sw = 'spectral_files/sp_sw_30_sub_nep_final_mstar'

os.system('mkdir '+directory)
S_start = 65.
S_inc = 5.
S_vec = []
i = 0

S_value = S_start
start = True
while True:

    S_value += S_inc

    if S_value < 114 and start:
        i = i+1
        S_vec.append(S_value)
        continue
    else:
        start = False
    fname = f'olr_{i}.nc'
    input_file = f'olr_{i-1}.nc'

    nml = f90nml.read('input.nml')

    nml['io_nml']['output_file'] = directory + '/' + fname
    
    if os.path.isfile(directory + '/' + input_file):
        nml['io_nml']['init_from_file'] = True
        nml['io_nml']['input_file'] = directory + '/' + input_file
    else:
        nml['io_nml']['init_from_file'] = False

    nml['radiation_nml']['Finc'] = S_value

    nml['surface_nml']['sensible_heat'] = False
    nml['timestep_nml']['accelerate'] = True
    
    if i==0:
        nml['convection_nml']['conv_switch'] = False
    else:
        nml['convection_nml']['conv_switch'] = True

    nml.write('input.nml', force=True)

    nml = f90nml.read('input_soc.nml')

    nml['socrates_rad_nml']['stellar_constant'] = S_value
    nml['socrates_rad_nml']['lw_spectral_filename'] = spectral_file_lw
    nml['socrates_rad_nml']['sw_spectral_filename'] = spectral_file_sw

    nml.write('input_soc.nml', force=True)

    command = ['python', 'run.py'] 

    proc = subprocess.run(command)

    if proc.returncode==99:
        while True:
            
            nml = f90nml.read('input.nml')
            tstep = nml['timestep_nml']['del_time']

            if tstep/2 > 8640:
                print('Reducing timestep to aid convergence')
                nml['timestep_nml']['del_time'] = tstep/2
                nml.write('input.nml', force=True)
            else:
                print('still stuck even with reduced timestep')
                sys.exit()
            proc = subprocess.run(command)
            if proc.returncode != 99:
                break

    if proc.returncode==100:
        while True:

            # Reduce increment in S

            if S_inc > 1.01:
                print('HERE', S_value, S_inc)
                S_value = S_value - np.ceil(S_inc/2)
                print('New S_value', S_value)
                S_inc = np.floor(S_inc/2)
                nml = f90nml.read('input.nml')
                nml['radiation_nml']['Finc'] = S_value
                nml.write('input.nml', force=True)

                nml = f90nml.read('input_soc.nml')
                nml['socrates_rad_nml']['stellar_constant'] = S_value
                nml.write('input_soc.nml', force=True)
                
            else:
                print('Runaway limit reached')
                sys.exit()

            proc = subprocess.run(command)
                
            if proc.returncode != 100:
                break
            
    if proc.returncode==101:
        while True:
            nml = f90nml.read('input.nml')
            nml['timestep_nml']['accelerate'] = False
            nml['surface_nml']['sensible_heat'] = True

            nml.write('input.nml', force=True)

            proc = subprocess.run(command)
            print('Return code: ', proc.returncode)
            if proc.returncode==0:
                nml = f90nml.read('input.nml')
                nml['timestep_nml']['accelerate'] = True
                nml['surface_nml']['sensible_heat'] = False
                nml['timestep_nml']['del_time'] = 86400
                nml.write('input.nml', force=True)
                break
            else:
                # Try with smaller timestep
                nml = f90nml.read('input.nml')
                if nml['timestep_nml']['del_time'] > 8640.:
                    nml['timestep_nml']['del_time'] = nml['timestep_nml']['del_time']/2
                    nml.write('input.nml', force=True)
                else:
                    print('Still stuck with reduced timestep')
                    sys.exit()
                    
                #test = xr.open_dataset(directory+'/'+fname+'.nc')
    #toa = (test.fup[0] + test.s_up[0] - test.fdn[0] - test.s_dn[0])/S
    #ts = test.Ts.data

    
    #if (np.absolute(toa) > 0.1):
    #    test.close()
    #    os.system('rm '+directory+'/'+fname+'.nc')
    #    i_end = i-1
    #    break
    #i_end = i
    S_vec.append(S_value)
    i+=1

S_vector = np.array(S_vec)
np.savetxt(directory+'/S_vec.txt', S_vector)

os.system('python merge_data.py '+ directory + ' 60 61')
