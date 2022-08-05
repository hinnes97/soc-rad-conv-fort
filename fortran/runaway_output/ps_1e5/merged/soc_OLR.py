import os
import nctools
import xarray as xr
import params
import subprocess
import numpy as np

def calc_SW_soc(p, T, q_h2o):
    """Do SW calculation with SOCRATES to test how much downwelling radiation
       penetrates to the radiative layer"""

    # Need values at mid-points of layers
    p_lay = (p[1:]+ p[:-1])/2
    t_lay = np.interp(p_lay, p, T)
    q_lay = np.interp(p_lay, p, q_h2o)
    # Mass fractions of all the species
    
    q_dry = 1 - q_lay
    q_ch4 = params.ch4_frac*q_dry
    q_h2he = q_dry - q_ch4
    q_h2 = 0.9*q_h2he
    q_he = 0.1*q_h2he

    tsurf = T[-1] # hope this is right
    p_surf = p[-1]


    surface_albedo=0.0
    solar_zenith_angle=0.0
    solar_toa = 1368./4.
    
    nctools.ncout_surf('profile.surf', 0,0,1,surface_albedo)
    nctools.ncout2d('profile.tstar',0,0,tsurf,'tstar',longname="Surface Temperature",units='K')
    nctools.ncout2d('profile.pstar',0,0,p_surf,'pstar',longname="Surface Pressure",units='PA')
    nctools.ncout3d('profile.t',0,0,p_lay,t_lay,'t',longname="Temperature",units='K')
    nctools.ncout3d('profile.tl',0,0,p,T,'tl',longname="Temperature",units='K')
    nctools.ncout3d('profile.p',0,0,p_lay,p_lay,'p',longname="Pressure",units='PA')
    nctools.ncout3d('profile.q',0,0,p_lay,q_lay,'q',longname="q",units='PPMM')
#    nctools.ncout3d('profile.h2o',0,0,pres_list,q_mr_list,'h2o',longname="h2o",units='PPMM')
    nctools.ncout3d('profile.h2',0,0,p_lay,q_h2,'H2',longname="Hydrogen",units='PPMM') # Add
    nctools.ncout3d('profile.he',0,0,p_lay,q_he,'He',longname="Helium",units='PPMM') # Add
    nctools.ncout3d('profile.ch4',0,0,p_lay,q_ch4,'CH4',longname="Methane",units='PPMM') # Add
    nctools.ncout2d('profile.szen',0,0,solar_zenith_angle,'szen',longname="Solar zenith angle",units='Degrees')
    nctools.ncout2d('profile.stoa',0,0,solar_toa,'stoa',longname="Solar Irradiance at TOA",units='WM-2')

    basename="profile"
    spec = params.spec_sw
    # Run Socrates:
    command = ("Cl_run_cdf -B", basename, "-s", spec, "-R 1 30 -ch 30 -S -g 8 0 -u -C 5 -t 15 -r")
    command = " ".join(command)

    os.system(command)
    os.system("fmove "+basename+ " currentsw")
    # We only actually care about OLR (also maybe contribution function)
    # Note last two coordinates are size 1 (holders for latitude and longitude)
    ds_netsw  =xr.open_dataset('currentsw.nflx')
    ds_swup = xr.open_dataset('currentsw.uflx')
    sw_net = np.sum(ds_netsw['nflx'].data, axis=0)[:,0,0]
    sw_up = np.sum(ds_swup['uflx'].data, axis=0)[:,0,0]
    print('SWUP', np.amax(sw_up))
    
#    ds_contrib = xr.open_dataset('currentlw.cff')
#    contrib = np.sum(ds_contrib['cff'].data, axis=0)[:,0,0]

    ds_netsw.close()
    ds_swup.close()
#    ds_contrib.close()

    return sw_net

def calc_OLR_soc(p, T, q_lay, p_lay, t_lay):
    """Calculate the OLR with SOCRATES radiation code"""

    # Need values at mid-points of layers
    #p_lay = (p[1:]- p[:-1])/np.log(p[1:]/p[:-1])
    #t_lay = np.interp(p_lay, p, T)
    #q_lay = np.interp(p_lay, p, q_h2o)
    # Mass fractions of all the species
    
    q_dry = 1 - q_lay
    q_ch4 = params.ch4_frac*q_dry
    q_h2he = q_dry - q_ch4
    q_h2 = 0.74812*q_h2he
    q_he = 0.25188*q_h2he

    tsurf = 436.7 # hope this is right
    p_surf = p[-1]

    
    surface_albedo=0.0
    nctools.ncout_surf('profile.surf', 0,0,1,surface_albedo)
    nctools.ncout2d('profile.tstar',0,0,tsurf,'tstar',longname="Surface Temperature",units='K')
    nctools.ncout2d('profile.pstar',0,0,p_surf,'pstar',longname="Surface Pressure",units='PA')
    nctools.ncout3d('profile.t',0,0,p_lay,t_lay,'t',longname="Temperature",units='K')
    nctools.ncout3d('profile.tl',0,0,p,T,'tl',longname="Temperature",units='K')
    nctools.ncout3d('profile.p',0,0,p_lay,p_lay,'p',longname="Pressure",units='PA')
    nctools.ncout3d('profile.pl',0,0,p,p, 'pl', longname='Pressure Levels', units='PA')
    nctools.ncout3d('profile.q',0,0,p_lay,q_lay,'q',longname="q",units='PPMM')
#    nctools.ncout3d('profile.h2o',0,0,pres_list,q_mr_list,'h2o',longname="h2o",units='PPMM')
    nctools.ncout3d('profile.h2',0,0,p_lay,q_h2,'H2',longname="Hydrogen",units='PPMM') # Add
    nctools.ncout3d('profile.he',0,0,p_lay,q_he,'He',longname="Helium",units='PPMM') # Add
    nctools.ncout3d('profile.ch4',0,0,p_lay,q_ch4,'CH4',longname="Methane",units='PPMM') # Add 
    basename="profile"
    spec = params.spec
    print(spec)
    # Run Socrates
    command = ("Cl_run_cdf -B", basename, "-s", spec, "-R 1 30 -ch 30 -I -g 4 -u -C 5 -t 12 -q -z 2 -x params.txt")
    command = " ".join(command)

    os.system(command)
    os.system("fmove "+basename+ " currentlw")
    # We only actually care about OLR (also maybe contribution function)
    # Note last two coordinates are size 1 (holders for latitude and longitude)
    ds_uflx  =xr.open_dataset('currentlw.uflx')
    OLR = np.sum(ds_uflx['uflx'].data, axis=0)[0,0,0]
    

    #ds_dflx = xr.open_dataset('currentlw.dflx')
    #test=  np.sum(ds_dflx['dflx'].data, axis=0)[:,0,0]
    
    #new = test[q_h2o>params.Rc*T/params.L/params.pi]

    #sb = 5.67e-8
    #bb = sb*T[q_h2o>params.Rc*T/params.L/params.pi]**4
    #print(np.sum(ds_uflx['uflx'].data, axis=0)[-1,0,0], test[-1], sb*T[-1]**4)

    #print('BB comparison: ', np.amax(np.r_[0,np.absolute((new-bb)/bb)]), (np.r_[0,new.T])[-1], (np.r_[0,bb.T])[-1], new-bb)
#    ds_contrib = xr.open_dataset('currentlw.cff')
#    contrib = np.sum(ds_contrib['cff'].data, axis=0)[:,0,0]

    ds_uflx.close()
#    ds_contrib.close()

    return OLR, 0
