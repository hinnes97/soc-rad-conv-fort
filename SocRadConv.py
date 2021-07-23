'''
MDH 28/01/19
Socrates radiative-convective model
'''

from timeit import default_timer as timer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import signal
matplotlib.use('Agg')
import SocRadModel
from atmosphere_column import atmos
import math,phys
import moistad

# # Font settings
# import matplotlib.pylab as pylab
# import matplotlib.font_manager as fm
# font = fm.FontProperties(family = 'Helvetica', fname = '/Users/tim/Dropbox/work/matplotlib_fonts/Helvetica/Helvetica.ttf')

"""------------------------------------------------------------------------ """
"""------------------------Thermodynamic constants------------------------- """
"""------------------------------------------------------------------------ """
"""
R = phys.water.R                         # J/(kg*K) specific gas constant of water vapor
Rcp = phys.water.Rcp                     # cp in J/(kg*K) specific heat constant of water vapor
L = phys.water.L_vaporization            # J/kg, latent heat of condensation of water vapor at 300K
esat = phys.satvps_function(phys.water)  # Saturation vapor pressure, arguments=T,T0,e0,MolecularWeight,LatentHeat
Tref = 350.                              # Reference temperature
pref = esat(Tref)                        # Reference pressure
"""

R = 461.9173324615943                 # J/(kg*K) specific gas constant of water vapor
Rcp = 0.2500905968931209              # cp in J/(kg*K) specific heat constant of water vapor
L = 2493000.0                         # J/kg, latent heat of condensation of water vapor at 300K
Tref = 350.                           # Reference temperature
pref = 46794.80349022859              # Reference pressure

""" Moist adjustment switch """
Moist_Adjustment = True

#Dew point temperature
def Tdew(p):
    return Tref/(1-(Tref*R/L)*math.log(p/pref)) 

def cond(T,p):

    eps = phys.H2O.MolecularWeight/phys.H2.MolecularWeight
    
    Lvap_R = phys.H2O.L_vaporization_TriplePoint*phys.H2O.MolecularWeight/phys.Rstar
    Lsub_R = phys.H2O.L_sublimation*phys.H2O.MolecularWeight/phys.Rstar
    Tref = phys.H2O.TriplePointT
    pref = phys.H2O.TriplePointP

    qsat = lambda L_R : pref*np.exp(-L_R*(1/T - 1/Tref))/p * eps

    if T>Tref:
        qsats = qsat(Lvap_R)
    else:
        qsats = qsat(Lsub_R)

    return qsats
    #q = np.minimum(q0, qsats)

    # Cold trap ensures q above cold point doesn't exceed value at coldest point
    #index = signal.argrelmin(T,order=2)[0][0]
    #q[:index] = q[index]
    
    #return q
    
def surf_Planck_nu(atm):
    h = 6.63e-34
    c = 3.0e8
    kb = 1.38e-23
    B = np.zeros(len(atm.band_centres))
    c1 = 1.191042e-5
    c2 = 1.4387752
    for i in range(len(atm.band_centres)):
        nu = atm.band_centres[i]
        B[i] = (c1*nu**3 / (np.exp(c2*nu/atm.ts)-1))

    B = B * atm.band_widths/1000.0
    return B

def planck_wn(wav, T):
        h = 6.626e-34
        c = 3.0e+8
        k = 1.38e-23
        a = 2e8*h*c*c*wav**3
        b = 100.0*h*c*wav/(k*T)
        intensity = a/ (  (np.exp(b) - 1.0) )
        return intensity

def RadConvEqm(output_dir, time_current, n_steps, Ts, stellar_toa_heating, p_s, h2o_ratio, co2_ratio, h2_ratio, ch4_ratio, co_ratio, n2_ratio, o2_ratio, he_ratio, n2o_ratio):
    #--------------------Set radmodel options-------------------
    #---Instantiate the radiation model---

    atm = atmos()

    #---Set up pressure array (a global)----
    atm.ps = p_s
    pstart = .995*atm.ps
    rat = (atm.ptop/pstart)**(1./atm.nlev)
    logLevels = [pstart*rat**i for i in range(atm.nlev+1)]
    logLevels.reverse()
    levels = [atm.ptop + i*(pstart-atm.ptop)/(atm.nlev-1) for i in range(atm.nlev+1)]
    atm.pl = np.array(logLevels)
    atm.p = (atm.pl[1:] + atm.pl[:-1]) / 2
    atm.dp = atm.pl[1:] - atm.pl[:-1]
    atm.q0 = h2o_ratio

    #==============Now do the calculation====================================

    atm.ts = 450#Ts
    
    atm.Rcp = phys.H2.Rcp
    atm.temp = atm.ts*(atm.p/atm.p[-1])**atm.Rcp  #Initialize on an adiabat
    atm.temp  = np.where(atm.temp<230.,230.,atm.temp)

    atm.temp = np.loadtxt('output/Ts7/atm_TP_profile1990.dat')[:,1]
    mintemp = signal.argrelmin(atm.temp, order=20)[0][-1] # Highest pressure min in temp
    atm.temp[:mintemp] = atm.temp[mintemp] # Sets min T to cold point

    atm.ts = atm.temp[-1]
    #atm.ts = np.loadtxt('../exp2/output/Ts0/atm_dt.dat')[-1,1]
    # atm.n_species = 2
    atm.n_species = 7

    Moist_adiabat=[Tdew(pp) for pp in atm.p*100.] # Tdew(p) with p in Pa
    #print("Moist_adiabat = " + str(Moist_adiabat))
    # # Water vapour
    # atm.mixing_ratios[0] = 1.e-5
    # # CO2
    # atm.mixing_ratios[1] = 1.e-5

    atm.mixing_ratios[0] = h2o_ratio # H2O
    atm.mixing_ratios[0] = np.loadtxt('output/Ts7/atm_h2o_profile1990.dat')[:,1]#h2o_ratio*np.ones(len(atm.p))
    atm.mixing_ratios[1] = co2_ratio # CO2
    atm.mixing_ratios[2] = 1-atm.mixing_ratios[0]  # H2
    atm.mixing_ratios[3] = ch4_ratio # CH4
    atm.mixing_ratios[4] = co_ratio  # CO
    atm.mixing_ratios[5] = n2_ratio  # N2
    atm.mixing_ratios[6] = o2_ratio  # O2
    atm.mixing_ratios[7] = n2o_ratio # N2O

    # Initialise previous OLR and TOA heating to zero
    PrevOLR = 0.
    PrevMaxHeat = 0.
    PrevTemp = 0.*atm.temp[:]

    # Arrays to keep track of Ts and the OLR in the loop
    #n_steps = 100 ! redefined in the running script
    Ts_array = np.ones(n_steps)
    OLR_array = np.ones(n_steps)
    dt_array = np.zeros(n_steps)

    #---------------------------------------------------------
    #--------------Initializations Done-----------------------
    #--------------Now do the time stepping-------------------
    #---------------------------------------------------------
    matplotlib.rc('axes',edgecolor='k')
    counter=0
    for i in range(n_steps):

        counter+=1
        atm = steps(atm, stellar_toa_heating)
        print("--- %s seconds ---" % (timer() - start_time)) # execution time

        if i % 5 == 0:
            print("Iteration", i, end =", ")
            #if 1==2:
                #plt.figure(figsize=(7,4))
                #plt.semilogy(atm.temp,atm.p)
                #plt.gca().invert_yaxis()
                #plt.ylabel('Pressure [mb]')
                #plt.xlabel('Temperature [K]')
                #plt.gca().xaxis.label.set_color('white')
                #plt.tick_params(axis='x', colors='white')
                #plt.gca().yaxis.label.set_color('white')
                #plt.tick_params(axis='y', colors='white')
                # plt.show()
            #print("OLR " + str(atm.LW_flux_up[-1]))
            print("OLR change " + str(atm.LW_flux_up[0]-PrevOLR))
            print("Max heating " + str(np.max(atm.total_heating)))
            print("Max dT " + str(abs(np.max(atm.temp-PrevTemp[:]))))

            print("OLR = " + str(PrevOLR)+" W/m^2,", "Max heating = " + str(np.max(atm.total_heating)))
            print("Skin T = ", atm.temp[0])

        # Reduce timestep if heating not converging
        #if abs(np.max(atm.temp-PrevTemp[:])) > 0.05 or abs(atm.temp[0]-atm.temp[1]) > 3.0:
        #    print("reducing timestep")
        #    atm.dt  = atm.dt*0.99
        #elif atm.dt < 1.0:
        #    print("increasing timestep")
        #    atm.dt = atm.dt/0.99

            
        #if abs(atm.LW_flux_up[0]-PrevOLR) < 0.1 and abs(np.max(atm.temp-PrevTemp[:])) < 0.5:
            #print("break")
            #print(PrevTemp[:]-atm.temp)
            #break    # break here

        PrevOLR = atm.LW_flux_up[0]
        PrevMaxHeat = abs(np.max(atm.total_heating))
        PrevTemp[:] = atm.temp[:]


        # Store the Ts and OLR values for each timestep
        Ts_array[i]  = atm.ts
        OLR_array[i] = PrevOLR
        dt_array[i]=atm.dt

        if i % 10 ==0:
            out_a = np.column_stack( ( atm.p, atm.temp, Moist_adiabat ) )
            np.savetxt( output_dir+f"atm_TP_profile{i}.dat", out_a )

            out_a = np.column_stack( ( atm.net_flux[:-1], atm.flux_down[:-1], atm.flux_up[:-1], atm.SW_flux[:-1], atm.LW_flux[:-1], atm.LW_flux_up[:-1], atm.contrib_LW_flux_up ) )
            np.savetxt( output_dir+f"atm_flux_profile{i}.dat", out_a )

            out_a = np.column_stack( ( atm.total_heating, atm.sw_heating, atm.lw_heating ) )
            np.savetxt( output_dir+f"atm_heating_profile{i}.dat", out_a )

            out_a = np.column_stack( ( atm.band_centres, atm.LW_spectral_flux_up[:,0]/atm.band_widths ) )
            np.savetxt( output_dir+f"atm_LWup_flux{i}.dat", out_a )

            out_a = np.column_stack( ( atm.band_centres, atm.LW_spectral_flux_net[:,0]/atm.band_widths ) )
            np.savetxt( output_dir+f"atm_LWnet_flux{i}.dat", out_a )

            out_a = np.column_stack( ( atm.band_centres, atm.SW_spectral_flux_net[:,0]/atm.band_widths ) )
            np.savetxt( output_dir+f"atm_SWnet_flux{i}.dat", out_a )

            out_a = np.column_stack( (atm.p, atm.mixing_ratios[0]) )
            np.savetxt( output_dir+f"atm_h2o_profile{i}.dat", out_a)
        
    out_a = np.column_stack( ( dt_array, Ts_array ) )
    np.savetxt( output_dir+"atm_dt.dat", out_a )
    
    return

# Dry adjustment routine
def dryAdj(atm):
    T = atm.temp
    p = atm.p
    dp = atm.dp
    #Rcp is a global
    #Downward pass
    for i in range(len(T)-1):
        T1,p1,dp1 = T[i],p[i],dp[i]
        T2,p2,dp2 = T[i+1],p[i+1],dp[i+1]
        pfact = (p1/p2)**atm.Rcp
        if T1 < T2*pfact:
            Tbar = (dp1*T1+dp2*T2)/(dp1+dp2) #Equal layer masses
                              #Not quite compatible with how
                              #heating is computed from flux
                              # FIXED 
            T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
            T1 = T2*pfact
            atm.temp[i] = T1
            atm.temp[i+1] = T2
    #Upward pass
    for i in range(len(T)-2,-1,-1):
        T1,p1,dp1 = T[i],p[i],dp[i]
        T2,p2,dp2 = T[i+1],p[i+1],dp[i+1]
        pfact = (p1/p2)**atm.Rcp
        if T1 < T2*pfact:
            Tbar = (dp1*T1+dp2*T2)/(dp1+dp2) #Equal layer masses
                              #Not quite compatible with how
                              #heating is computed from flux
            T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
            T1 = T2*pfact
            atm.temp[i] = T1
            atm.temp[i+1] = T2

#Moist adjustment routine.
def moistAdj(atm):
    # HII changed to make this dilute moist adiabat adjustment
    T = atm.temp
    p = 100*atm.p # Need to unit convert to Pa from mbar
    dp = 100*atm.dp
    q0 = atm.q0

    #Leconte 2018 composition effect
    vap = phys.H2O
    vap.update()
    omega = 1 - phys.H2.MolecularWeight/vap.MolecularWeight

    # Cold trap
    qsats = np.array([cond(T,100*p) for (p,T) in zip(atm.p,atm.temp)])
    mins = signal.argrelmin(qsats,order=20)[0]
    if mins.size>0:
        indices = mins[q0>qsats[mins]]
        if indices.size>0:
            index=indices[-1]
            atm.mixing_ratios[0][:index+1] = qsats[index]
        else:
            index = 0
    else:
        index = 0
            
    for i in np.arange(index,len(T)-1):
        qsat = cond(atm.temp[i+1],100*atm.p[i+1]) #Unit conversion into Pa

        if q0 > qsat:
            atm.mixing_ratios[0][i+1] = qsat
            T1,p1,dp1 = T[i],p[i],dp[i]
            T2,p2,dp2 = T[i+1],p[i+1],dp[i+1]
            pfact = (p1/p2)**moistad.dlntdlnp(np.log(p2), np.log(T2),phys.H2,phys.H2O)

            if T2 > vap.TriplePointT:
                L = vap.L_vaporization
            else:
                L = vap.L_sublimation

            qcrit = vap.R*T2/L/omega
            
            if T1 < T2*pfact and atm.mixing_ratios[0][i+1] < qcrit:
                Tbar = (dp1*T1+dp2*T2)/(dp1+dp2) #Equal layer masses
                #Not quite compatible with how
                #heating is computed from flux
                # FIXED 
                T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
                T1 = T2*pfact
                atm.temp[i] = T1
                atm.temp[i+1] = T2
                atm.mixing_ratios[0][i+1] = cond(T2,p2)
                atm.mixing_ratios[0][i]   = cond(T1,p1)
        else:
            atm.mixing_ratios[0][i+1] = q0
                
                
    #Upward pass
    for i in np.arange(len(T)-2,index-1,-1):
        qsat = cond(atm.temp[i+1],100*atm.p[i+1])
        if q0>qsat:
            atm.mixing_ratios[0][i+1] = qsat
            T1,p1,dp1 = T[i],p[i],dp[i]
            T2,p2,dp2 = T[i+1],p[i+1],dp[i+1]
            pfact = (p1/p2)**moistad.dlntdlnp(np.log(p2), np.log(T2),phys.H2,phys.H2O)
            if T2 > vap.TriplePointT:
                L = vap.L_vaporization
            else:
                L = vap.L_sublimation

            qcrit = vap.R*T2/L/omega

            if T1 < T2*pfact and atm.mixing_ratios[0][i+1] < qcrit:
                Tbar = (dp1*T1+dp2*T2)/(dp1+dp2) #Equal layer masses
                #Not quite compatible with how
                #heating is computed from flux
                T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
                T1 = T2*pfact
                atm.temp[i] = T1
                atm.temp[i+1] = T2
                atm.mixing_ratios[0][i+1] = cond(T2,p2)
                atm.mixing_ratios[0][i]   = cond(T1,p1)
        else:
            atm.mixing_ratios[0][i+1] = q0
    
#Define function to do time integration for n steps
def steps(atm, stellar_toa_heating):
    # HII remove excess H2O here
    #atm.mixing_ratios[0] = cond(atm.mixing_ratios[0],atm.temp,atm.p, atm.q0)
    
    atm = SocRadModel.radCompSoc(atm, stellar_toa_heating)
    dT = atm.total_heating*atm.dt
    #Limit the temperature change per step
    dT = np.where(dT>5.,5.,dT)
    dT = np.where(dT<-5.,-5.,dT)
    #Midpoint method time stepping
    #changed call to r.  Also modified to hold Tg fixed
    atm = SocRadModel.radCompSoc(atm, stellar_toa_heating)
    dT = atm.total_heating*atm.dt
    #Limit the temperature change per step
    dT = np.where(dT>5.,5.,dT)
    dT = np.where(dT<-5.,-5.,dT)
    atm.temp += dT
    #
    dTmax = max(abs(dT)) #To keep track of convergence

    #   Do the surface balance
    pturb = 1000. # Depth of boundary layer near surface, over which turbulent heating will be spread. In mb
             #Note if this is made too large, it can suppress convection
    iturb = np.searchsorted(atm.p,atm.p[-1]-pturb)
    dpTurb = atm.p[-1]-atm.p[iturb]
    pturbT = atm.p[iturb] #Top of boundary layer, mb
    kturb = 20. # W/m**2/K
    F_geo = 0. # Geothermal flux [W/m2] (up to 1 W/m2)
    mu = 1e6 # Thermal inertia of the ground [J.m−2.K−1]
    dp = atm.p[-1]-atm.p[-2]
    g = 15.58
    mua = phys.H2.cp*100.*dpTurb/g #Low layer atmospheric thermal inertia
    da = 24.*3600.
    #atm.temp[atm.p>pturbT] += -atm.dt*da*kturb*(atm.temp[-1] - atm.ts)/mua


    #atm.ts += (atm.dt*da/mu)*(stellar_toa_heating - atm.net_flux[-1] - kturb*(atm.ts - atm.temp[-1]) + F_geo )
    #atm.ts += (atm.dt*da/mu)*(-atm.net_flux[-1] - kturb*(atm.ts - atm.temp[-1]) + F_geo )

    # HII conserve energy at the bottom boundary
    fluxup = atm.flux_down[-1]
    atm.ts = (fluxup/phys.sigma)**0.25

    # HII Trying to replicate flux boundary condition
    print(atm.net_flux[-1])
    #Dry adjustment step
    for iadj in range(10):
        dryAdj(atm)
        moistAdj(atm)

    # After all the H2O change need to adjust H2 before next SOC call
    atm.mixing_ratios[2] = 1 - atm.mixing_ratios[0]

    Tad = atm.temp[-1]*(atm.p/atm.p[-1])**atm.Rcp
    #** Temporary kludge to keep stratosphere from getting too cold
    atm.temp = np.where(atm.temp<50.,50.,atm.temp)  #**KLUDGE
    #
    #Dummies for separate LW and stellar. **FIX THIS**
    fluxStellar = fluxLW = heatStellar = heatLW = np.zeros(atm.nlev)
    return atm
