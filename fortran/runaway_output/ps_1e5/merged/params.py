import phys

N = 1000
pup = 0.1
#kappa =0.01
#kappa=0.01
kappa = 0.007742840658352238 # Found by fitting to the exact Socrates value of OLR in steam limit
k2 = 0.00 # Binary coefficient for h2 CIA, need to fit this value to a SOCRATES OLR curve
kap0 = 1.6e-4#2.5e-4
grav = 9.8
tau_0 = 0.0
dry = phys.H2
wet = phys.H2O
L = wet.L_vaporization
Rc = wet.R
pi = 1 - dry.MolecularWeight/wet.MolecularWeight
mass_path = 1.e4
strat_temp= 0.
CI = True
pure_radiative = False
# SW heating
kappa_sw = 1.e-5 # kg/m^2, works roughly


# SOCRATES
ch4_frac=0.0# Proportion of dry air that is methane (by mass)
#spec='spec_files/sp_lw_300_sub_nep'
spec='spec_files/sp_lw_300_sub_nep'
spec_sw = 'spec_files/sp_sw_300_sub_nep_sol'
soc = False

#
old_scheme = False
qmin = 1
rad_strat= False
