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

plt.semilogy(ds.Tf.data[conv_mask], ds.pfull.data[conv_mask], lw=2, color='C0')
plt.semilogy(ds.Tf.data[~conv_mask], ds.pfull.data[~conv_mask], color='C0')
plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')

plt.gca().invert_yaxis()
twin=plt.gca().twiny()
twin.loglog(ds.q.data[0,:], ds.pfull.data, color='C1')
twin.set_xlabel('Specific Humidity')

plt.tight_layout()

plt.savefig(args.output)
