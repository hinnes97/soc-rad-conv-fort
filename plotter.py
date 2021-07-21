import numpy as np
import matplotlib.pyplot as plt

sig = 5.67e-8
#for i in range(29,30):
with open('data/state_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    pe = data[:,0]
    fup = data[:,2]
    fd = data[:,3]
    sd = data[:,4]
with open('data/state_FINAL.csv') as fh:
    data = np.genfromtxt(fh)
    p = data[:,0]
    T = data[:,1]
        
def analytic(p, taulwinf,tauswinf):
    tau = taulwinf*(p/p[-1])
    gamma = tauswinf/taulwinf
    S0 = 1368/4

    # To test dryadj put Fint=100
    Fint=0
    
    test = Fint*(tau+1)+ S0*(1+1/gamma + (gamma-1/gamma)*np.exp(-gamma*tau))
    return (test/2/sig)**0.25

def F_p(tau, gamma, S0):
    return S0/2*(1+1/gamma + np.exp(-gamma*tau)*(1-1/gamma))

def F_m(tau,gamma,S0):
    return F_p(tau,gamma,S0) - S0*np.exp(-gamma*tau)

plt.figure()

tau = 256*pe/pe[-1]
gamma = 0.04
S0=1368/4

plt.semilogy(F_p(tau,gamma,S0), pe, F_m(tau,gamma,S0), pe, S0*np.exp(-gamma*tau), pe,  F_p(tau,gamma,S0)-F_m(tau,gamma,S0)-S0*np.exp(-gamma*tau), pe)

plt.semilogy(fup, pe, fd, pe, sd, pe, fup-fd-sd,pe, ls = '--')
plt.gca().invert_yaxis()

plt.figure()
plt.semilogy(T, p)
plt.semilogy(analytic(p,256, 256*0.04), p)
plt.semilogy(T[-1]*(p/p[-1])**(2/7), p)
plt.gca().invert_yaxis()

plt.show()
