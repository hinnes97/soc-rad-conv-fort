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
import moistad as ma
import convection as conv
import phys


class atmos:
    def __init__(self, p_top, p_s, Ne, T_init,S0):
        self.pe = np.logspace(np.log10(p_top), np.log10(p_s), Ne)
        self.Tf = T_init
        self.Ne = Ne
        self.Nf = Ne - 1
        dp = np.diff(self.pe)
        self.dp = np.r_[dp[0], dp]

        self.S0 = S0
        self.Fint =100 #Use 100 to test convective adjustment
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
        
        # Convection
        self.N0=None
        self.Rcp = 2/7
        self.dry = phys.H2
        self.wet = phys.H2O

        self.dry_mask = np.full(len(self.T), False)
        self.wet_mask = np.full(len(self.T), False)

        #Moisture
        self.q0 = 0.1
        self.qsat = ma.satq(self.T, self.p, self.dry, self.wet)
        self.q = np.minimum(self.q0, self.qsat)
        self.q = conv.cold_trap(self.T, self.p, self.q)[1]
        
        # Initialise residual
        self.R = self.calc_residual(self.Te)
        
    def interp_to_edge(self,Tf, pf, pe):

        Te = np.zeros(self.Ne)

        logpf = np.log(pf)
        logpe = np.log(pe)
       
        f = spi.CubicSpline(logpf, Tf, extrapolate=True)

        return f(logpe)

    def calc_residual(self, T):
        self.f_down = rad.ir_flux_down(T,self.pe)

        # Need to calculate all f_up even if N0!=None, because need upwards flux boundary
        # condition from bottom of the convective layer

        self.f_up = rad.ir_flux_up(T,self.pe)
        
        self.s_down = rad.sw_flux_down(self.S0, self.pe)

        return  self.s_down + self.f_down - self.f_up + self.Fint

    def calc_jacobian(self):
        dT = 1 # Lower this for more accuracy?

        # Calculate regions we actually want to do radiation on (exclude adiabats), but
        # include beginning of adiabatic regions so that region can grow/shrink
        
        dry_blocks = conv.start_stop(self.dry_mask, True)
        wet_blocks = conv.start_stop(self.wet_mask, True)

        new_dry_mask = np.copy(self.dry_mask)
        new_wet_mask = np.copy(self.wet_mask)
        
        for m,n in dry_blocks:
            new_dry_mask[m] = False
        for m,n in wet_blocks:
            new_wet_mask[m] = False

        rad_mask = ~(new_dry_mask | new_wet_mask)
        len_rad = np.sum(rad_mask)

        jacob = np.zeros((len_rad, len_rad))
        dummy_q = np.zeros_like(self.q)
        
        for i,j in enumerate(np.arange(self.Ne)[rad_mask]):#range(len(jacob[0])):
            T_dash = np.copy(self.T)
            T_dash[j] += dT

            # Currently deciding against changing temperatures below convective zone - reasoning
            # that when temperature finally incremented we don't know for sure that the convective
            # adjustment will cause the levels below it to also change (i.e. convective layer could
            # shrink -- see if it works
            
            Te_dash = self.interp_to_edge(T_dash, self.p, self.pe)

            Fdash_i = self.calc_residual(Te_dash)[rad_mask]
            
            jacob[:,i] = (Fdash_i - self.R[rad_mask])/dT

        self.jacob = jacob
        self.rad_mask = rad_mask

    def update_state(self):
        max_dT = 5
        self.calc_jacobian()
        self.dT = np.linalg.solve(self.jacob, -self.R[self.rad_mask])
        
        # Limit magnitude of dT
        self.dT[self.dT>5] = max_dT
        self.dT[self.dT<-5] = -max_dT

        self.T[self.rad_mask] += self.dT

        self.T[self.T<150] = 150
        self.Tf  = self.T[1:]

        self.Te = self.interp_to_edge(self.T, self.p, self.pe)

        # Calculate residual
        self.R[self.rad_mask] = self.calc_residual(self.Te)[self.rad_mask]

        print(f'Residual = {self.R[self.rad_mask]}')
        #print(f'Max residual = {np.amax(np.absolute(self.R[self.rad_mask])):.2e} W/m^2')
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

    def dump_state(self, path, stamp):

        # Save state at one timestamp into one file
        with open(path+"_"+stamp+".csv", 'wb') as fh:
            np.savetxt(fh, np.transpose([self.pf, self.Tf]))

        with open(path+'_fluxes_'+stamp+'.csv', 'wb') as fh:
            np.savetxt(fh, np.transpose([self.pe, self.Te, self.f_up, self.f_down, self.s_down,
                                         self.q, self.dry_mask, self.wet_mask]))
                        
    def run_to_equilib(self, m, n, path):
        for j in range(n):
            print(f'n = {j} -------------------------------')
            for i in range(m):
                #self.calc_residual(self.Te)
                self.dump_state(path, str(i))
                self.update_state()
                if np.amax(np.absolute(self.R[:self.N0]))<1e-10:
                    break
            
            self.T[1:-1] = 0.25*self.T[:-2] + 0.5*self.T[1:-1] + 0.25*self.T[2:]
            self.T[0] = 0.75*self.T[0]+0.25*self.T[1]
            #self.T[-1] = 0.75*self.T[-1] + 0.25*self.T[-2]

            # Perform moist and dry adjustment
            #self.T, self.dry_mask = conv.dry_adjust(self.T, self.p, self.dry_mask)
            #self.T, self.q, self.wet_mask = conv.moist_adjust(self.T, self.p, self.dp, self.q,
            #                                                    self.wet_mask, whole_atm = True, n_iter=4)

            self.Te = self.interp_to_edge(self.T, self.p, self.pe)
            # If moist adjustment overwrites dry adjustment, adjust dry mask
            self.dry_mask[self.wet_mask] = False
            
        self.dump_state(path, 'FINAL')
        
        
if __name__ == '__main__':
    pt = 1e2
    ps = 1e6
    Ne = 100


    
    pp = np.logspace(1,5,Ne-1)
    def analytic(p, taulwinf,tauswinf):
        tau = taulwinf*(p/p[-1])
        gamma = tauswinf/taulwinf
        S0 = 1368/16

        test = S0*(1+1/gamma + (gamma-1/gamma)*np.exp(-gamma*tau))
        return (test/2/rad.sig)**0.25

    #t_init = analytic(pp, 10,6)
    t_init=np.linspace(200,400,Ne-1)
    atm = atmos(pt, ps, Ne, t_init, 1368/16)

    atm.run_to_equilib(100, 5, 'data/state')
