import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline


sig = 5.670374419e-8
taulwinf = 10
tauswinf = 6

def flux_down_integrand(tau_dash, tau, spline):
    return sig*spline(np.log(tau_dash))**4*np.exp(-(tau-tau_dash))

def flux_up_integrand(tau_dash, tau, spline):
    return sig*spline(np.log(tau_dash))**4*np.exp(-(tau_dash-tau))

def ir_flux_down(Te, pe):
    tau = ir_tau(pe, taulwinf)

    spline = CubicSpline(np.log(tau), Te, extrapolate=True)
    
    I_m = np.zeros_like(Te)
    for i in range(1, len(I_m)):
        I_m[i] = quad(flux_down_integrand, tau[0], tau[i], args=(tau[i], spline))[0]

    return I_m

def ir_flux_up(Te,pe):
    tau = ir_tau(pe, taulwinf)

    spline = CubicSpline(np.log(tau), Te, extrapolate=True)
    
    I_p = np.zeros_like(Te)
    I_p[-1] = sig*Te[-1]**4
    
    for i in range(len(I_p)-2,-1,-1):
        I_p[i] = I_p[i+1]*np.exp(-(tau[-1] - tau[i])) + quad(flux_up_integrand, tau[i], tau[-1], args=(tau[i], spline))[0]

    return I_p

def sw_flux_down(S0, p):
    tau_sw = sw_tau(p, tauswinf)
    return S0*np.exp(-tau_sw)

def ir_tau(p, tauinf):
    return tauinf*(p/p[-1])

def sw_tau(p, tauinf):
    return tauinf*(p/p[-1])
    
