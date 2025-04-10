import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import xarray as xr
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile1',type=str)
parser.add_argument('infile2',type=str)
parser.add_argument('outfile', type=str)
args = parser.parse_args()

ds1 = xr.open_dataset(args.infile1)
ds2 = xr.open_dataset(args.infile2)
plt.figure()

for i,ds in enumerate([ds1,ds2]):
    conv_mask = ds.dry_mask.data>0

    plt.semilogy(ds.Tf.data[conv_mask], ds.pfull.data[conv_mask], lw=2, color='C'+str(i))
    plt.semilogy(ds.Tf.data[~conv_mask], ds.pfull.data[~conv_mask], color='C'+str(i))
    
plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')

plt.gca().invert_yaxis()
twin=plt.gca().twiny()

for i,ds in enumerate([ds1,ds2]):
    twin.loglog(ds.q.data[0,:], ds.pfull.data, color='C'+str(i), ls='--')
twin.set_xlabel('Specific Humidity')

plt.tight_layout()

plt.savefig(args.outfile)
