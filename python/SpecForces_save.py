import numpy as np
import glob as glob
import matplotlib.cm as cm
import Analysis as aa
from scipy.io import netcdf

##########################
### Control parameters ###
##########################
# Choose input directory
idir = '../'
# Choose output directory
odir = '../force_test_mode1/'
# Name for  saved file
xavg = False
incomp = True
if xavg:
    name = 'spec_xavgpost.bak'
else:
    name= 'spec.bak'

# Spectral truncation mode numbers
nx_c = 16
nz_c = 64
# Run for which to calculate this
run = 'Lx16Lz64_Re66_IC66'
# Number of saved restress file
n = -1 # Taking the last one saved
# Reynolds stress type
Retype = '2d'

######################################################################
# First let's load the Reynolds stresses to calculate the force
restress = sorted(glob.glob(idir+run+'/restress_'+Retype+'*'))

if incomp:
    # Extract the forces in spectral space:
    time, Re, Lx, Lz, specF_full = aa.SpecReForces(restress[n],Retype=Retype)
    # Project onto incompressible field:
    specF_comp,specF = aa.div_comp(Lx,Lz,specF_full)
else:
    # Extract the forces in spectral space:
    time, Re, Lx, Lz, specF = aa.SpecReForces(restress[n],Retype=Retype)
    
# Coarse grain by spectral truncation
N_new = nz_c
MM_new = nx_c
M_new = 2*(MM_new-1)
KK = specF.shape[2]
specF_filt = np.zeros((N_new,M_new,KK),dtype=complex)
specF_filt[:,:MM_new,:] = -specF[:N_new,:MM_new,:] # Minus sign because now we're going to feed in <u'.grad(u')>, as a nonlinear term, not a forcing on the RHS
specF_filt[:,MM_new:,:] = -specF[:N_new,-MM_new+2:,:]

# Set modes to zero:
specF_filt[:,:,0] = 0.0
# specF_filt[:,:,1] = 0.0
specF_filt[:,:,2] = 0.0
specF_filt[:,:,3] = 0.0
# specF_filt[:,:,4] = 0.0
specF_filt[:,:,5] = 0.0
specF_filt[:,:,6] = 0.0
specF_filt[:,:,7] = 0.0
# specF_filt[:,:,8] = 0.0
specF_filt[:,:,9] = 0.0
specF_filt[:,:,10] = 0.0

if xavg:
    # Set m=/=0 to zero.
    specF_filt[:,1:,:] = 0.0 + 1.j*0.0

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
    