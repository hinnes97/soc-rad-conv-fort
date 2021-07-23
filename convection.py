import numpy as np
import phys
import moistad as ma
from scipy import signal

dry = phys.H2
wet = phys.H2O
q0 = 0.1
Rcp = 2./7.

def start_stop(a, trigger_val):
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False,np.equal(a, trigger_val),False]

    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])

    # Get the start and end indices with slicing along the shifting ones
    return [x for x in zip(idx[::2], idx[1::2])]

def dry_adjust(T, p, dry_mask):

    dlnT = np.log(T)
    dlnp = np.log(p)

    grad = np.gradient(dlnT, dlnp)

    # Add small number to dry mask so that code notes regions already on adiabat
    new_mask = np.copy(dry_mask)
    dry_blocks = start_stop(dry_mask, True)

    for m,n in dry_blocks:
        new_mask[m] = False

    grad[new_mask] += 0.01

    for k in range(len(dlnT)-1):
        if grad[k] > Rcp:
            dry_mask[k:k+2] = True
            T[k+1] = T[k]*(p[k+1]/p[k])**Rcp
        else:
            dry_mask[k:k+2] = False
#            if np.all(grad[k:-1]> self.Rcp):
#                self.N0=k+1
#                self.T[k:] = self.T[k]*(self.p[k:]/self.p[k])**self.Rcp
#                break
    return T, dry_mask

def moist_adjust(Tin, p, dp, q, wet_mask, whole_atm=True, n_iter=1):
    T = np.copy(Tin)

    # If adjusting entire atmosphere, calculate cold trap. If not, we are doing as part of the
    # radiative calculation, so set cold trap index to 0
    if whole_atm == True:
        cti, q = cold_trap(T,p,q)
    else:
        cti = 0

    #grad = np.gradient(dlnT, dlnp)
    #if self.N0 != None:
    #    grad[self.N0:] += 0.001

    # Need to do downwards/upwards pass unlike dry adjustment because gradient changes with
    # changing temperature

    for n in range(n_iter):
        # Downwards pass
        for k in range(cti, len(T) - 1):
            qsat1 = ma.satq(T[k], p[k], dry, wet)
            qsat2 = ma.satq(T[k+1], p[k+1], dry, wet)

            if q0>qsat1 and q0 > qsat2:
                q[k] = qsat1
                q[k+1] = qsat2

                T1,p1,dp1 = T[k], p[k],dp[k]
                T2,p2,dp2 = T[k+1], p[k+1], dp[k+1]
                pfact = (p1/p2)**ma.dlntdlnp(np.log(p2), np.log(T2),dry,wet)

                if T2 > wet.TriplePointT:
                    L = wet.L_vaporization
                else:
                    L = wet.L_sublimation

                new_mask = np.copy(wet_mask)
                moist_blocks = start_stop(wet_mask, True)

                for m,n in moist_blocks:
                    new_mask[m] = False

                if new_mask[k] == True:
                    test = 1-0.01
                else:
                    test = 1

                print(T2*pfact/T1)
                if test < T2*pfact/T1:
                    Tbar = (dp1*T1+dp2*T2)/(dp1+dp2) #Equal layer masses                    
                    T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
                    T1 = T2*pfact
                    T[k] = T1
                    T[k+1] = T2
                    q[k+1] = ma.satq(T2,p2, dry, wet)
                    q[k]   = ma.satq(T1,p1, dry, wet)

                    wet_mask[k:k+2] = True
                else:
                    wet_mask[k:k+2] = False

            else:
                if qsat1 > q0:
                    q[k] = q0
                if qsat2 > q0:
                    q[k+1] = q0

                wet_mask[k:k+2] = False

        # Upwards pass
        for k in range(len(T) - 2, cti-1, -1):
            qsat1 = ma.satq(T[k], p[k], dry, wet)
            qsat2 = ma.satq(T[k+1], p[k+1], dry, wet)

            if q0>qsat1 and q0>qsat2:
                q[k] = qsat1
                q[k+1] = qsat2

                T1,p1,dp1 = T[k], p[k], dp[k]
                T2,p2,dp2 = T[k+1], p[k+1], dp[k+1]
                pfact = (p1/p2)**ma.dlntdlnp(np.log(p2), np.log(T2),dry,wet)

                if T2 > wet.TriplePointT:
                    L = wet.L_vaporization
                else:
                    L = wet.L_sublimation

                new_mask = np.copy(wet_mask)
                moist_blocks = start_stop(wet_mask, True)

                for m,n in moist_blocks:
                    new_mask[m] = False

                if new_mask[k] == True:
                    test = 1-0.01
                else:
                    test = 1

                print(T2*pfact/T1)
                if test < T2*pfact/T1:
                    Tbar = (dp1*T1+dp2*T2)/(dp1+dp2) #Equal layer masses                    
                    T2 = (dp1+dp2)*Tbar/(dp2+dp1*pfact)
                    T1 = T2*pfact
                    T[k] = T1
                    T[k+1] = T2
                    q[k+1] = ma.satq(T2,p2,dry, wet)
                    q[k]   = ma.satq(T1,p1, dry,wet)
                    wet_mask[k:k+2] = True
                else:
                    wet_mask[k:k+2] = False
            else:
                if qsat1 > q0:
                    q[k] = q0
                if qsat2 > q0:
                    q[k+1] = q0
                wet_mask[k:k+2] = False            

    return T, q, wet_mask

def cold_trap(T, p, q):
    qsat = ma.satq(T,p,dry,wet)

    mins = signal.argrelmin(qsat,order=20)[0]
    if mins.size>0:
        indices = mins[q0>qsat[mins]]
        if indices.size>0:
            index=indices[-1]
            q[:index+1] = qsat[index]
        else:
            index = 0
    else:
        index = 0

    return index, q
