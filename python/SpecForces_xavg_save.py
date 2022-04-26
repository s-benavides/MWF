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
odir = '../force_test_5mode/' # Taking modes off based on forcing.
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
specF_filt[:,0,:] = -specF[:N_new,:] # Minus sign because now we're going to feed in <u'.grad(u')>, as a nonlinear term, not a forcing on the RHS

# Setting modes to zero:
# specF_filt[:,:,0] = 0.0 # v1, 9 modes
# specF_filt[:,:,4] = 0.0

# specF_filt[:,:,2] = 0.0 # v2, 4 modes
# specF_filt[:,:,3] = 0.0 
# specF_filt[:,:,5] = 0.0
# specF_filt[:,:,6] = 0.0
# specF_filt[:,:,7] = 0.0 
# specF_filt[:,:,8] = 0.0
# specF_filt[:,:,10] = 0.0

specF_filt[:,:,0] = 0.0 # v3, 5 modes
specF_filt[:,:,3] = 0.0 
specF_filt[:,:,4] = 0.0
specF_filt[:,:,6] = 0.0
specF_filt[:,:,9] = 0.0 
specF_filt[:,:,10] = 0.0

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