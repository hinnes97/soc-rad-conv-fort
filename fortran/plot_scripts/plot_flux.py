import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import xarray as xr
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile',type=str)
parser.add_argument('outfile', type=str)
args = parser.parse_args()

ds = xr.open_dataset(args.infile)
plt.figure()
conv_mask = ds.dry_mask.data>0

fig, ax = plt.subplots(1,2, layout='constrained', sharey=True)
conv_data = np.where(conv_mask, ds.Tf.data, np.nan)
ax[0].semilogy(ds.Tf.data[conv_mask], ds.pfull.data[conv_mask], lw=5, color='C0')
ax[0].semilogy(ds.Tf.data, ds.pfull.data, color='C0')
ax[0].set_xlabel('Temperature [K]')
ax[0].set_ylabel('Pressure [Pa]')

ax[0].invert_yaxis()
twin=ax[0].twiny()
twin.loglog(ds.q[0,:].data, ds.pfull.data, color='C1')
twin.set_xlabel('Specific Humidity')

# Fluxes
ax[1].plot(-ds.s_dn.data, ds.pedge.data, label='s_dn')
ax[1].plot(ds.s_up.data, ds.pedge.data, label='s_up')
ax[1].plot(-ds.fdn.data, ds.pedge.data, label='f_dn')
ax[1].plot(ds.fup.data, ds.pedge.data, label='f_up')
ax[1].plot(ds.turb_flux.data, ds.pedge.data, label='turb')
ax[1].plot(ds.fup.data+ds.s_up.data - ds.fdn.data - ds.s_dn.data + ds.turb_flux.data, ds.pedge.data , label = 'net')
ax[1].legend()

ax[1].set_xscale("symlog", linthresh=0.01)
ax[1].set_xlabel('Flux (W/m^2)')

plt.savefig(args.outfile)
