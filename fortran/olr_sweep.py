import numpy as np
import os
import xarray as xr
import subprocess
import f90nml
import sys

S_vec = np.arange(50.,150.,5.)
S_vec = np.r_[np.arange(50.,130.,5.)]

directory = 'final_mstar_new_new'
spectral_file_lw = 'spectral_files/sp_lw_30_sub_nep_final'
spectral_file_sw = 'spectral_files/sp_sw_30_sub_nep_final_mstar'

os.system('mkdir '+directory)
S_start = 50.
S_inc = 5.
S_vec = []
i = 0
tstep_init = 86400.
S_value = S_start
start = True
tstep = tstep_init
pure_timestep = False

counter = 0
while True:
#    if S_value < 114 and start:
#        i = i+1
#        S_vec.append(S_value)
#        S_value += S_inc
#        continue
    #else:
    #    start = False
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

    if pure_timestep:
        nml['timestep_nml']['accelerate'] = False
        nml['surface_nml']['sensible_heat'] = True
    else:
        nml['timestep_nml']['accelerate'] = True
        nml['surface_nml']['sensible_heat'] = False
    
    if i==0:
        nml['convection_nml']['conv_switch'] = False
        tstep = tstep_init
    else:
        nml['convection_nml']['conv_switch'] = True

    nml['timestep_nml']['del_time'] = tstep
    
    nml.write('input.nml', force=True)

    nml = f90nml.read('input_soc.nml')

    nml['socrates_rad_nml']['stellar_constant'] = S_value
    nml['socrates_rad_nml']['lw_spectral_filename'] = spectral_file_lw
    nml['socrates_rad_nml']['sw_spectral_filename'] = spectral_file_sw

    nml.write('input_soc.nml', force=True)

    command = ['python', 'run.py'] 

    proc = subprocess.run(command)

    if proc.returncode==99 or proc.returncode ==100:
        nml = f90nml.read('input.nml')
        tstep = nml['timestep_nml']['del_time']

        if tstep > 8640:
            print('Reducing timestep to aid convergence')
            tstep = tstep/2
            continue
        else:
            if S_inc > 0.26:
                
                S_inc_new = S_inc/2

                if S_inc_new < 0.24:
                    S_inc_new = 0.25
                    
                S_value = S_value - S_inc + S_inc_new

                S_inc = S_inc_new
                tstep = tstep_init
                continue
            else:
                print('Runaway limit hit')
                sys.exit()

    if proc.returncode==101:
        counter +=1
        print( 'COUNTER = ', counter)
        if counter>2:
            pure_timestep = True
            print('PURE TIMESTEPPING ON')
            continue
        nml = f90nml.read('input.nml')
        nml['timestep_nml']['accelerate'] = False
        nml['surface_nml']['sensible_heat'] = True

        nml.write('input.nml', force=True)

        tstep_old = tstep
        while tstep > 8640.:
            proc = subprocess.run(command)
            print('Return code: ', proc.returncode)
        
            if proc.returncode==0:
                nml = f90nml.read('input.nml')
                nml['timestep_nml']['accelerate'] = True
                nml['surface_nml']['sensible_heat'] = False
                tstep = tstep_old
                nml.write('input.nml', force=True)
                break
            else:
                tstep = tstep/2
                nml = f90nml.read('input.nml')
                nml['timestep_nml']['del_time'] = tstep
                nml.write('input.nml', force=True)
        else:
            print('PURE TIMESTEPPING SHOULD BE ON')
            pure_timestep = True
            continue
            
    S_vec.append(S_value)
    i+=1
    S_value += S_inc
    tstep = tstep_init

    if S_value>150:
        break
S_vector = np.array(S_vec)
np.savetxt(directory+'/S_vec.txt', S_vector)

os.system('python merge_data.py '+ directory + ' 60 61')

def reduce_increment(S_value, S_inc):
    S_inc_new = S_inc/2

    S_value = S_value - S_inc + S_inc_new
    S_inc = S_inc_new

    return S_value, S_inc
