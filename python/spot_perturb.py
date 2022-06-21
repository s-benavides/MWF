import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import matplotlib.cm as cm
import Analysis as aa
from scipy.io import netcdf
import scipy.ndimage as ndim

# Random seed:
rng = np.random.default_rng(12345)

# Amplitude of perturbation:
amp = 1.

# Choose input directory
idir = '../'

# Name of original file to use
spot_IC = idir+'spot_test_Re70_512_amp_10_outs/state0011.cdf.dat'

# Load file
time,Re,Lx,Lz,mp = aa.read_cdf(spot_IC,intype= 'mpt')

for ii in range(1,11):
    # Choose output directory
    odir = '../spot_test_Re70_512_ens_'+str(ii)+'/'

    # Name for  saved file
    name = 'IC.bak'
#     name = 'state0001.cdf.dat'
    
    # Perturb
    mpt = np.copy(mp)*(1+amp*(rng.random(size=mp.shape)-0.5)) # so that random num. between 1-0.5*amp and 1+0.5*amp

    ##############
    ### Saving ###
    ##############
    print("Saving...")
    with netcdf.NetCDFFile(odir+name,'w') as f:
        # Create dimensions
        f.createDimension('N', 512) #NN
        f.createDimension('M', 1022) #M
        f.createDimension('K',7) #K
        f.createDimension('ReIm',2)

        # define attributes
        f.t = 0.0
        f.t_step = 0
        f.Re = 70.0
        f.alpha = 2*np.pi/Lx
        f.gamma = 2*np.pi/Lz

        # Define variable
        mpt_in = f.createVariable('mpt', 'd', ('ReIm','N','M','K'))
        mpt_in.K=7
        mpt_in.MM=512
        mpt_in.NN=512
        mpt_in[:] = mpt

    print("Done! Saved %s in %s" % (name,odir))
