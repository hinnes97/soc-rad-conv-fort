"""
Class intended to solve the Rad-Conv problem given fluxes calculated in another module. Trying 
to keep it as general as possible such that it will still work if grey gas replaced with some
other scheme (e.g. picket-fence, Socrates)
"""

# General imports
#--------------------
import numpy as np
import scipy.interpolate as spi
import os
import matplotlib.pyplot as plt

# Below import your specific radiation/convection scheme
#--------------------------------------------------------
#import grey as rad
import short_char as rad

class atmos:
    def __init__(self, p_top, p_s, Ne, T_init,S0):
        self.pe = np.logspace(np.log10(p_top), np.log10(p_s), Ne)
        self.Tf = T_init
        self.Ne = Ne
        self.Nf = Ne - 1

        self.S0 = S0
        self.Fint =0# rad.sig*80**4
        self.f_up=np.ones_like(self.pe)
        self.f_down=np.ones_like(self.pe)
        self.s_down = np.ones_like(self.pe)
        
        # Define pf levels similaraly to FMS
        self.pf = (self.pe[1:] - self.pe[:-1])/(np.log(self.pe[1:]) - np.log(self.pe[:-1]))
        self.Te = self.interp_to_edge(self.Tf, self.pf, self.pe)

        # Have to include temp at TOA in temp vector to have equal number of fluxes + temps
        self.T = np.concatenate((self.Te[0], self.Tf), axis=None)
        self.p = np.concatenate((self.pe[0], self.pf), axis=None)
        self.dT = np.zeros_like(self.Tf)
        self.R = np.zeros_like(self.T)

        self.N0=None

    def interp_to_edge(self,Tf, pf, pe):

        Te = np.zeros(self.Ne)

        logpf = np.log(pf)
        logpe = np.log(pe)
        
        f = spi.CubicSpline(logpf, Tf, extrapolate=True)#, fill_value="extrapolate")

        return f(logpe)

    def calc_residual(self, T):
        self.f_down[:self.N0] = rad.ir_flux_down(T,self.pe[:self.N0])
        self.f_up[:self.N0] = rad.ir_flux_up(T,self.pe[:self.N0])
        self.s_down[:self.N0] = rad.sw_flux_down(self.S0, self.pe[:self.N0])

        return  self.s_down[:self.N0] + self.f_down[:self.N0] - self.f_up[:self.N0] + self.Fint

    def calc_jacobian(self):
        dT = 1 # Lower this for more accuracy?

        if not self.N0:
            jacob = np.zeros((self.Ne, self.Ne))
        else:
            jacob = np.zeros((self.N0, self.N0))

        # Calculate residual
        self.R[:self.N0] = self.calc_residual(self.Te[:self.N0]) 
        
        for j in range(len(jacob[0])):
            T_dash = np.copy(self.T)

            T_dash[j] += dT

            Te_dash = self.interp_to_edge(T_dash, self.p, self.pe)

            Fdash_j = self.calc_residual(Te_dash[:self.N0])
            
            jacob[:,j] = (Fdash_j - self.R[self.N0])/dT

        self.jacob = jacob

    def update_state(self):
        max_dT = 5
        self.calc_jacobian()
        self.dT = np.linalg.solve(self.jacob, -self.R[:self.N0])

        # Smooth dT
        #self.dT[1:-1] = 0.25*self.dT[:-2] + 0.5*self.dT[1:-1] + 0.25*self.dT[2:]
        
        # Limit magnitude of dT
        self.dT[self.dT>5] = max_dT
        self.dT[self.dT<-5] = -max_dT
        print(f'Max residual = {np.amax(np.absolute(self.R)):.2e} W/m^2')

        self.T[:self.N0] += self.dT

        self.T[self.T<150] = 150
        self.Tf  = self.T[1:]

        self.Te = self.interp_to_edge(self.T, self.p, self.pe)

        #### Below for timestepping version
        #max_dT = 5
        #const = 0.1
        #self.R = self.calc_residual(self.Te)
        #
        #self.dT = -(self.R[1:] - self.R[:-1])*const
        #
        #self.dT[self.dT>max_dT] = max_dT
        #self.dT[self.dT<-max_dT] = -max_dT
        #
        #self.R = self.calc_residual(self.interp_to_edge(self.Tf+0.5*self.dT,self.pf, self.pe))
        #
        #self.dT = -(self.R[1:] - self.R[:-1])*const
        #        
        #self.dT[self.dT>max_dT] = max_dT
        #self.dT[self.dT<-max_dT] = -max_dT
        #
        #self.Tf += self.dT
        #self.Te = self.interp_to_edge(self.Tf, self.pf, self.pe)
        #
        #print(np.amax(np.absolute(self.dT)), np.amax(np.absolute(self.R)))

    # Note convective adjustment is NOT working at the moment
    def dry_adjust(self):
        Rcp=2/7
        dlnT = np.log(self.T)
        dlnp = np.log(self.p)

        grad = np.gradient(dlnT, dlnp)
        if self.N0 != None:
            grad[self.N0+1:] += 0.001

        self.N0=None
        for k in range(self.Ne):
            if np.all(grad[k:-1]> Rcp):
                self.N0=k
                break

        self.T[self.N0:] = self.T[self.N0]*(self.p[self.N0:]/self.p[self.N0])**Rcp
        self.Te = self.interp_to_edge(self.T, self.p, self.pe)
                                
    def dump_state(self, path, stamp):

        # Save state at one timestamp into one file
        with open(path+"_"+stamp+".csv", 'wb') as fh:
            np.savetxt(fh, np.transpose([self.pf, self.Tf]))

        with open(path+'_fluxes_'+stamp+'.csv', 'wb') as fh:
            np.savetxt(fh, np.transpose([self.pe, self.Te, self.f_up, self.f_down, self.s_down]))
                        
    def run_to_equilib(self, m, n, path):
        for j in range(n):
            for i in range(m):
                #self.calc_residual(self.Te)
                self.dump_state(path, str(i))
                self.update_state()
                if np.amax(np.absolute(self.R[:self.N0]))<1e-10:
                    break
            
            self.T[1:-1] = 0.25*self.T[:-2] + 0.5*self.T[1:-1] + 0.25*self.T[2:]
            self.T[0] = 0.75*self.T[0]+0.25*self.T[1]
            #self.T[-1] = 0.75*self.T[-1] + 0.25*self.T[-2]

            #self.dry_adjust()
            #print(self.N0)
            #print(self.R[:self.N0])
        
        self.dump_state(path, 'FINAL_noadj')
        
        
if __name__ == '__main__':
    pt = 1e2
    ps = 1e6
    Ne = 100

    pp = np.logspace(1,5,Ne-1)
    def analytic(p, taulwinf,tauswinf):
        tau = taulwinf*(p/p[-1])
        gamma = tauswinf/taulwinf
        S0 = 1368/4

        test = S0*(1+1/gamma + (gamma-1/gamma)*np.exp(-gamma*tau))
        return (test/2/rad.sig)**0.25

    #t_init = analytic(pp, 10,6)
    t_init=np.linspace(200,400,99)
    atm = atmos(pt, ps, Ne, t_init, 1368/4)

    atm.run_to_equilib(100, 5, 'data/state')
