import numpy as np

sig = 5.670374419e-8
taulwinf = 10
tauswinf = 6

def ir_flux_down(Te, pe):
    
    tau = ir_tau(pe, taulwinf)
    dtau = np.diff(tau)

    # Define constants for interpolation of source
    lc = linear_coeff(dtau)
    pc = parab_coeff(dtau)

    I_m = np.zeros_like(Te)
    S = sig*Te**4

    # Parabolic for every term but last (linear)
    dI_m = np.zeros((len(I_m) - 1))
    #dI_m = lc.a_m*S[:-1] + lc.b_m*S[1:]
    
    dI_m[:-1] = pc.a_m*S[:-2] + pc.b_m*S[1:-1] + pc.g_m*S[2:]
    dI_m[-1] = lc.a_m[-1]*S[-2] + lc.b_m[-1]*S[-1]
    
    for k in range(1,len(I_m)):
        I_m[k] = I_m[k-1]*np.exp(-dtau[k-1]) + dI_m[k-1]
        
    return I_m

def ir_flux_up(Te, pe):
    
    tau = ir_tau(pe, taulwinf)
    dtau = np.diff(tau)

    lc = linear_coeff(dtau)
    pc = parab_coeff(dtau)

    S = sig*Te**4
    
    I_p = np.zeros_like(Te)
    I_p[-1] = S[-1]

    dI_p = np.zeros((len(I_p) - 1))
    #dI_p = lc.b_p*S[:-1] + lc.g_p*S[1:]
    
    dI_p[1:] = pc.a_p*S[:-2]+ pc.b_p*S[1:-1] + pc.g_p*S[2:]
    dI_p[0] = lc.b_p[0]*S[0] + lc.g_p[0]*S[1]

    for k in range(len(I_p)-2,-1,-1):
        I_p[k] = I_p[k+1]*np.exp(-dtau[k]) + dI_p[k]
        
    return I_p

def sw_flux_down(S0, p):
    tau_sw = sw_tau(p, tauswinf)
    return S0*np.exp(-tau_sw)

def ir_tau(p, tauinf):
    return tauinf*(p/p[-1])

def sw_tau(p, tauinf):
    return tauinf*(p/p[-1])
    
class linear_coeff:
    def __init__(self, dtau):
        
        e0 = 1 - np.exp(-dtau)
        e1 = dtau - e0
        
        self.a_m = e0 - e1/dtau
        self.b_m = e1/dtau

        self.b_p = e1/dtau
        self.g_p = e0 - e1/dtau

class parab_coeff:
    def __init__(self, dtau):

        e0 = 1 - np.exp(-dtau)
        e1 = dtau - e0
        e2 = dtau**2 - 2*e1

        self.a_m = e0[:-1] + (e2[:-1] - (dtau[1:] + 2*dtau[:-1])*e1[:-1])/(dtau[:-1]*(dtau[1:] + dtau[:-1]))
        self.b_m = ((dtau[1:] + dtau[:-1])*e1[:-1] - e2[:-1])/(dtau[:-1]*dtau[1:])
        self.g_m = (e2[:-1] - dtau[:-1]*e1[:-1])/(dtau[1:]*(dtau[1:] + dtau[:-1]))

        self.a_p = (e2[1:] - dtau[1:]*e1[1:])/(dtau[:-1]*(dtau[1:] + dtau[:-1]))
        self.b_p = ((dtau[1:] + dtau[:-1])*e1[1:] - e2[1:])/(dtau[:-1]*dtau[1:])
        self.g_p = e0[1:] + (e2[1:] - (dtau[:-1] + 2*dtau[1:])*e1[1:])/(dtau[1:]*(dtau[1:] + dtau[:-1]))
