import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import matplotlib.cm as cm
import Analysis as aa
from scipy.io import netcdf
import scipy.ndimage as ndim

##########################
### Control parameters ###
##########################
# Choose input directory
idir = '../'
# Choose output directory
odir = '../force_test_xavgpre/'
# Name for  saved file
name = 'spec_xavgpre.bak'
# Spectral truncation mode numbers
nx_c = 16
nz_c = 64
# Run for which to calculate this
run = 'Lx16Lz64_Re66_IC66'
# Number of saved restress file
n = '2'

######################################################################
# Extract the forces in spectral space:
time,Re,Lx,Lz,z,specF = aa.SpecReForces_xavg(idir+run,n)
# Coarse grain by spectral truncation
N_new = nz_c
MM_new = nx_c
M_new = 2*(MM_new-1)
KK = int(specF.shape[1])
specF_filt = np.zeros((N_new,M_new,KK),dtype=complex)
specF_filt[:,0,:] = specF[:N_new,:]

with netcdf.NetCDFFile(odir+name,'w') as f:
    # Create dimensions
    f.createDimension('N', N_new)
    f.createDimension('M', M_new)
    f.createDimension('KK',KK)
    f.createDimension('ReIm',2)
    
    # define attributes
    f.t = 0.0
    f.t_step = 0
    f.Re = Re
    f.alpha = 2*np.pi/Lx
    f.gamma = 2*np.pi/Lz
    
    # Define variable
    force = f.createVariable('spec', 'd', ('ReIm','N','M','KK'))
    force.KK=KK
    force.MM=MM_new
    force.NN=N_new
    force[0,:,:,:] = np.real(specF_filt)
    force[1,:,:,:] = np.imag(specF_filt)
    
print("Done! Saved %s in %s" % (name,odir))