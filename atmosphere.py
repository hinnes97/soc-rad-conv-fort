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
#import short_char as rad
import grey as rad
import short_char as rad2

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

    def interp_to_edge(self,Tf, pf, pe):

        Te = np.zeros(self.Ne)

        logpf = np.log(pf)
        logpe = np.log(pe)
        
        f = spi.CubicSpline(logpf, Tf, extrapolate=True)#, fill_value="extrapolate")

        return f(logpe)

    def calc_residual(self, T):
        self.f_down = rad2.ir_flux_down(T,self.pe)
        self.f_up = rad2.ir_flux_up(T,self.pe)
        self.s_down = rad2.sw_flux_down(self.S0, self.pe)

        return  self.s_down + self.f_down - self.f_up + self.Fint

    def calc_jacobian(self):
        dT = 1 # Lower this for more accuracy?

        jacob = np.zeros((self.Ne, self.Ne))

        # Calculate residual
        self.R = self.calc_residual(self.Te) 
        
        for j in range(self.Ne):
            T_dash = np.copy(self.T)

            T_dash[j] += dT

            Te_dash = self.interp_to_edge(T_dash, self.p, self.pe)

            Fdash_j = self.calc_residual(Te_dash)
            
            jacob[:,j] = (Fdash_j - self.R)/dT

        self.jacob = jacob

    def update_state(self):
        max_dT = 5
        #self.calc_jacobian()
        #plt.imshow(np.log(np.absolute(self.jacob)))
        #plt.show()
        #self.dT = np.linalg.solve(self.jacob, -self.R)
        # Limit magnitude of dT
        #print(self.dT)
        const = 0.1
        self.R = self.calc_residual(self.Te)

        #plt.semilogy(self.f_up, self.pe)
        #plt.semilogy(self.f_down, self.pe)
        #plt.semilogy(self.s_down,self.pe)
        #plt.semilogy(self.R, self.pe)
        #plt.gca().invert_yaxis()
        #plt.show()
        
        self.dT = -(self.R[1:] - self.R[:-1])*const
        
        self.dT[self.dT>max_dT] = max_dT
        self.dT[self.dT<-max_dT] = -max_dT
        
        self.R = self.calc_residual(self.interp_to_edge(self.Tf+0.5*self.dT,self.pf, self.pe))

        self.dT = -(self.R[1:] - self.R[:-1])*const
                
        self.dT[self.dT>max_dT] = max_dT
        self.dT[self.dT<-max_dT] = -max_dT

        self.Tf += self.dT
        self.Te = self.interp_to_edge(self.Tf, self.pf, self.pe)

        print(np.amax(np.absolute(self.dT)), np.amax(np.absolute(self.R)))
        #self.Tf  = self.T[1:]

    def dump_state(self, path, stamp):

        # Save state at one timestamp into one file
        with open(path+"_"+stamp+".csv", 'wb') as fh:
            #np.savetxt(fh, np.transpose([self.p, self.T, self.dT, self.R, self.Te]))
            np.savetxt(fh, np.transpose([self.pf, self.Tf, self.dT]))

        with open(path+'_fluxes_'+stamp+'.csv', 'wb') as fh:
            np.savetxt(fh, np.transpose([self.pe, self.Te, self.f_up, self.f_down, self.s_down]))
                        
    def run_to_equilib(self, m, path):
        for i in range(m):
            #self.calc_residual(self.Te)
            self.dump_state(path, str(i))
            self.update_state()        

if __name__ == '__main__':
    pt = 1e1
    ps = 1e5
    Ne = 100

    pp = np.logspace(1,5,99)
    def analytic(p, taulwinf,tauswinf):
        tau = taulwinf*(p/p[-1])
        gamma = tauswinf/taulwinf
        S0 = 1368/4

        test = S0*(1+1/gamma + (gamma-1/gamma)*np.exp(-gamma*tau))
        return (test/2/rad.sig)**0.25

    #t_init = analytic(pp, 10,6)
    t_init=300*np.ones(99)   
    atm = atmos(pt, ps, Ne, t_init, 1368/4)

    atm.run_to_equilib(10000, 'data/state')
