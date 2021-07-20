import numpy as np
import matplotlib.pyplot as plt

sig = 5.67e-8
for i in range(1):
    with open('data/state_fluxes_'+str(i)+'.csv') as fh:
        data = np.genfromtxt(fh)
        print(data.shape)
        pe = data[:,0]
        fup = data[:,2]
        fd = data[:,3]
        sd = data[:,4]
    with open('data/state_'+str(i)+'.csv') as fh:
        data = np.genfromtxt(fh)
        p = data[:,0]
        T = data[:,1]
        


def analytic(p, taulwinf,tauswinf):
    tau = taulwinf*(p/p[-1])
    gamma = tauswinf/taulwinf
    S0 = 1368/4

    test = S0*(1+1/gamma + (gamma-1/gamma)*np.exp(-gamma*tau))
    return (test/2/sig)**0.25

def F_p(tau, gamma, S0):
    return S0/2*(1+1/gamma + np.exp(-gamma*tau)*(1-1/gamma))

def F_m(tau,gamma,S0):
    return F_p(tau,gamma,S0) - S0*np.exp(-gamma*tau)

plt.figure()

tau = 10*pe/pe[-1]
gamma = 6/10
S0=1368/4

plt.semilogy(F_p(tau,gamma,S0), pe, F_m(tau,gamma,S0), pe, S0*np.exp(-gamma*tau), pe,  F_p(tau,gamma,S0)-F_m(tau,gamma,S0)-S0*np.exp(-gamma*tau), pe)
plt.semilogy(fup, pe, fd, pe, sd, pe, fup-fd-sd,pe)
plt.gca().invert_yaxis()

plt.show()
