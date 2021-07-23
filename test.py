import numpy as np
import moistad as ma
import phys
import convection as conv
import matplotlib.pyplot as plt

def start_stop(a, trigger_val):
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False,np.equal(a, trigger_val),False]

    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])

    # Get the start and end indices with slicing along the shifting ones
    return [x for x in zip(idx[::2], idx[1::2]-1)]

#p = np.logspace(5,0,100)
#Ts = 250

#dry = 250*(p/p[0])**(2/7)

#moistad = np.zeros_like(dry)

#for i, (pp, tt) in enumerate(zip(p,dry)):
#    print('p,t', pp,tt)
#    moistad[i] = ma.dlntdlnp(np.log(pp), np.log(tt), phys.H2, phys.H2O)
#    print(moistad[i] - phys.H2.Rcp)

    
#tad,pad = ma.moistadiabat(p,Ts,phys.H2,phys.H2O)


#plt.loglog(tad,pad)
#plt.loglog(tad[0]*(pad/pad[0])**(phys.H2.Rcp), pad)
#plt.gca().invert_yaxis()
#plt.show()

pe = np.logspace(1,5,100)
pf = (pe[1:]-pe[:-1])/(np.log(pe[1:]) - np.log(pe[:-1]))
dp = np.diff(pe)

t = 280*(pf/pf[-1])**(3/7)
q = np.zeros_like(t)
wet_mask = np.full(len(t), False)

tnew, qnew, wet_mask = conv.moist_adjust(t,pf, dp, q,wet_mask,whole_atm=True)

print('NEW OUT', tnew)
print('OLD OUT', t)
plt.loglog(tnew, pf)
plt.loglog(t, pf)
plt.gca().invert_yaxis()
plt.show()
