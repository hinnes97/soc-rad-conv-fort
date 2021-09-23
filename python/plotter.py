import numpy as np
import matplotlib.pyplot as plt
import moistad as ma
import convection as conv
import phys

dry = phys.H2
wet = phys.H2O
sig = 5.67e-8
#for i in range(29,30):
with open('data/state_default_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    p = data[:,0]
    T = data[:,1]
    #fup = data[:,2]
    #fd = data[:,3]
    #sd = data[:,4]
    #q=data[:,5]
    #dry_mask=np.array([bool(x) for x in data[:,6]])
    #wet_mask=np.array([bool(x) for x in data[:,7]])

with open('data/state_increase_tau_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    Tir = data[:,1]

with open('data/state_decrease_tau_sw_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    Tsw = data[:,1]

with open('data/state_pbroad_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    Tbr = data[:,1]

with open('data/state_pbroad_high_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    Thigh = data[:,1]

with open('data/state_pbroad_Fint_high_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    Tfl = data[:,1]

with open('data/state_low_sw_tau_CIA_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    Tlswcia = data[:,1]

with open('data/state_twoSWband_fluxes_FINAL'+'.csv') as fh:
    data = np.genfromtxt(fh)
    T2sw = data[:,1]    


def analytic(p, taulwinf,tauswinf):
    tau = taulwinf*(p/p[-1])
    gamma = tauswinf/taulwinf
    S0 = 1368/4
    
    # To test dryadj put Fint=100
    Fint=5.67e-8*70**4
    print(taulwinf, tauswinf, S0, Fint)
    test = Fint*(tau+1)+ S0*(1+1/gamma + (gamma-1/gamma)*np.exp(-gamma*tau))
    return (test/2/sig)**0.25

def F_p(tau, gamma, S0):
    return S0/2*(1+1/gamma + np.exp(-gamma*tau)*(1-1/gamma))

def F_m(tau,gamma,S0):
    return F_p(tau,gamma,S0) - S0*np.exp(-gamma*tau)

plt.figure()

#tau = 256*p/p[-1]
#gamma = 0.04
#S0=1368/16

#plt.semilogy(F_p(tau,gamma,S0), p, F_m(tau,gamma,S0), p, S0*np.exp(-gamma*tau), p,  F_p(tau,gamma,S0)-F_m(tau,gamma,S0)-S0*np.exp(-gamma*tau), p)

#plt.semilogy(fup, p, fd, p, sd, p, fup-fd-sd,p, ls = '--')
#plt.gca().invert_yaxis()

#print(wet_mask)
#print(dry_mask)
def adiabat(Ts, p):
    return Ts*(p/p[-1])**phys.H2.Rcp

plt.figure()
plt.semilogy(T, p, label='Using this currently')
plt.semilogy(Tir, p, label='Infrared tau*4')
plt.semilogy(Tsw, p, label='SW tau/4')
plt.semilogy(Tbr, p, label='H2 CIA')
plt.semilogy(Thigh, p, label='H2 CIA high tau')
plt.semilogy(Tfl, p, label='H2 CIA high Fint (x2 K218b Tint)')
plt.semilogy(Tlswcia, p, label='H2 CIA low sw tau')
plt.semilogy(adiabat(Tfl[-1],p), p, label='Adiabat')
plt.semilogy(T2sw, p, label= 'Two SW bands')
plt.legend()
#plt.semilogy(analytic(p,4*256, 256*0.04), p)
#plt.semilogy(T[-1]*(p/p[-1])**(2/7), p)
#plt.scatter(T[wet_mask][0], p[wet_mask][0],marker='x',color='C0')
#plt.scatter(T[wet_mask][-1], p[wet_mask][-1], marker='x',color='C0')
#plt.scatter(T[dry_mask][0], p[dry_mask][0],marker='x',color='C1')
#plt.scatter(T[dry_mask][-1], p[dry_mask][-1], marker='x',color='C1')
#plt.plot(T[wet_mask][-1]*(p/p[wet_mask][-1])**(2/7), p)
plt.gca().invert_yaxis()

#tnew, new_dry_mask = conv.dry_adjust(T, p, np.r_[np.diff(p)[0], np.diff(p)],dry_mask, n_iter=1000)
#tnewer,q_new, new_wet_mask = conv.moist_adjust(T,p, np.r_[np.diff(p)[0], np.diff(p)], wet_mask,q, n_iter=1000)
#plt.loglog(tnew, p, ls='--')
#plt.loglog(tnewer, p, ls='--')

# Try and get moist adjustment to work on the layer that's convecting in the simulation
#tnew, qnew, new_wet_mask = conv.moist_adjust(T[wet_mask], p[wet_mask], np.r_[np.diff(p[wet_mask]),np.diff(p[wet_mask])[-1]], q[wet_mask], np.full(len(T[wet_mask]), True), whole_atm=True, n_iter=1000)
#plt.figure()
#plt.loglog(tnew, p[wet_mask])
#plt.loglog(T[wet_mask], p[wet_mask])
#plt.gca().invert_yaxis()
#print(wet_mask, new_wet_mask)
# True moist adiabat
#ttrue, ptrue = ma.moistadiabat(np.flip(p[wet_mask]), tnew[-1], phys.H2, phys.H2O)
#plt.loglog(ttrue, ptrue)

#qsat = ma.satq(T, p, dry,wet)
#plt.figure()
#plt.loglog(q, p, qwet, p)
#plt.loglog(q_new, p, qsat, p)
#plt.gca().invert_yaxis()
plt.show()
