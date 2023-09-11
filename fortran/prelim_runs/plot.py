import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

mpl.rcParams['text.usetex'] = True 
mpl.rcParams['text.latex.preamble'] = r'\usepackage[cm]{sfmath} \usepackage{siunitx} \usepackage{mhchem}'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'cm'


mpl.rcParams['font.size'] = 10
sns.set_palette('colorblind')
lw2 = 7.126625 # Exact double-column width in inches
lw = 3.365 # single column in inches
phi = 0.5*(5**0.5 + 1)
hgt = 9.129

shami_dat = np.genfromtxt('K2-18b-quiet-1kbar-1myr-steady-state.txt', skip_header=2)
with open ('K2-18b-quiet-1kbar-1myr-steady-state.txt') as f:
    l = f.readline()

    species = np.array(f.readline().split()[3:])
    
P = shami_dat[:,0]/10.0
T = shami_dat[:,1]

test = xr.open_dataset('fiducial.nc')


xs = shami_dat[:,3:]
mmw_list = np.array([17,2,18,1,16,13,12,14,15,16,24,26,25,27,28,29,30,28,44,31,30,29,31,32,43,32,42,41,14,15,26,27,30,16,28,17,30,29,31,32,31,28,43,46,44,50,30,29,31,44,47,42,33,34,51,41,40,53,39,38,40,51,53,78,77,16,16,14,4])
mmw = np.sum(xs*mmw_list, axis=1)

print(mmw)
qs = xs*mmw_list/mmw[:,np.newaxis]

filt = np.any(qs[:50] > 1.e-3, axis=0)
spec_main = species[filt]

inds = np.argsort(np.average(qs[:50,filt], axis=0))[::-1]
print('Most important species:')

for i,spec in enumerate(spec_main):

    print(f'{spec:<10}: {np.average(qs[:50,filt], axis=0)[inds][i]:7.2E}')

plt.figure(figsize=[lw, lw/1.2])
plt.semilogy(T,P, label='Helios, 100x Solar Metallicity')
plt.semilogy(test.Tf, test.pfull, label='Grey Gas')

plt.gca().invert_yaxis()
plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')
plt.title('K2-18 b, 1x Instellation')

plt.grid(alpha=0.5, ls='--')
plt.gca().set_xticks(np.arange(200.,1600,200.))
plt.tight_layout()
plt.savefig('TP_init.pdf')


q = 0.35
print('Average opacity:')
print(f'LW = {2.e-2*q + (1-q)*2.e-3:10.4E} m^2/kg')
print(f'SW = {4.e-4*q + (1-q)*2.e-5:10.4E} m^2/kg')

