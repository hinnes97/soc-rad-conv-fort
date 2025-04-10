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

plt.semilogy(ds.pt.data, ds.pfull.data, lw=2, color='C0', label='PT')
#plt.semilogy(ds.pt.data[~conv_mask], ds.pfull.data[~conv_mask], color='C0', label='PT')
plt.semilogy(ds.vpt.data, ds.pfull.data, lw=2, color='C1', label='VPT')
#plt.semilogy(ds.vpt.data[~conv_mask], ds.pfull.data[~conv_mask], color='C1', label='VPT')

plt.xlabel('Potential Temperature [K]')
plt.ylabel('Pressure [Pa]')

plt.gca().invert_yaxis()
plt.legend()
plt.tight_layout()

plt.savefig(args.outfile)
