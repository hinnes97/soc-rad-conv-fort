import numpy as np
import matplotlib.pyplot as plt
import phys
from scipy import integrate

def satvp(T,wet):
    """Calculates saturation vapour pressure"""
    T0 = wet.TriplePointT
    p0 = wet.TriplePointP

    # if T>wet.TriplePointT:
    #     L = wet.L_vaporization
    # else:
    #     L = wet.L_sublimation

    L = wet.L_vaporization
    return p0*np.exp( -L/wet.R*(1/T - 1/T0) )

def cold_trap(T, p, dry, wet, T_strat):
    
    
    arr = T<T_strat

    T[arr] = T_strat
    q = calc_q(T, p, dry, wet)
    try:
        index = np.argmin(q)
        q[:index] = q[index]
        #print(q[:index])
    except (ValueError, IndexError):
        q[arr] = q[arr][-1]
        
    #print('COLD TRAP', q[arr])
    return T, q
        
def calc_q(T, p,dry, wet):
    ratio = satvp(T, wet)/p
    eps = wet.MolecularWeight/dry.MolecularWeight
    
    return eps*ratio/(1 + (eps - 1)*ratio)

def dew_point(p, wet):
    """Finds the dew point temperature for a given pressure (inverts satvp)"""
    T0 = wet.TriplePointT
    p0 = wet.TriplePointP

    L = wet.L_vaporization

    return T0/(1 - wet.R*T0/L * np.log(p/p0))

def dlntdlnp(lnp,lnT,dry,wet):
    """Calculates moist adiabatic gradient (see Ding 2016)
       This value is correct, even in the nondilute regime
    """
    T = np.exp(lnT)
    p = np.exp(lnp)
    
    dry.update()
    wet.update()

    Rd,Rcpd,cpd = dry.R,dry.Rcp,dry.cp
    Rw,cpw = wet.R,wet.cp
    eps = wet.MolecularWeight/dry.MolecularWeight

    # if T>wet.TriplePointT:
    #     L = wet.L_vaporization
    # else:
    #     L = wet.L_sublimation

    L = wet.L_vaporization
    psat = satvp(T,wet)
    rsat = eps*psat/p
    pa = p - psat
    
    #rsat = qsat/(1-qsat)
    
    num = 1 + (L/Rd/T)*rsat
    den = 1 + ( (cpw/cpd) + ( (L/Rw/T) - 1 )*(L/cpd/T) )*rsat
    
    F = Rcpd*num/den

    dlnpdlnt = (psat/p)*L/Rw/T + pa/p/F
    
    return 1./dlnpdlnt

def moistadiabat(p,Ts,dry,wet):
    """Numerically integrates a moist adiabat from surface temp and pressure"""
    #ps = p[0]
    #pend = p[-1]
    
    ps = p[-1]
    pend = p[0]


    lnpend = np.log(pend)

    lnps,lnTs = np.log(ps),np.log(Ts)

    sol = integrate.solve_ivp(dlntdlnp,(lnps,lnpend),[lnTs],args=(dry,wet),t_eval=np.flip(np.log(p)), method='DOP853')
    #sol = integrate.solve_ivp(dlntdlnp,(lnps,lnpend),[lnTs],args=(dry,wet), t_)
    pad = np.flip(np.exp(sol.t))
    Tad = np.flip(np.exp(sol.y)[0])
    return Tad,pad

