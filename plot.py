import numpy as np
import matplotlib.pyplot as plt
import moistad as ma
import convection as conv
import phys

dry = phys.H2
wet = phys.H2O
sig = 5.67e-8

with open('data/test_FINAL.csv') as fh:
    data = np.genfromtxt(fh)
    p = data[:,0]
    T = data[:,1]

def adiabat(Ts, p):
    return Ts*(p/p[-1])**phys.H2.Rcp

plt.figure()
plt.semilogy(T, p, label='Using this currently')
plt.gca().invert_yaxis()
plt.show()
