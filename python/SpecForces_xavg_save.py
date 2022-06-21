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
odir = '../force_test_fstream_xavg/' # Taking modes off based on forcing.
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

##################################################
### x_avged force, and removing vertical modes ###
##################################################
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

# specF_filt[:,:,0] = 0.0 # v3, 5 modes
# specF_filt[:,:,3] = 0.0 
# specF_filt[:,:,4] = 0.0
# specF_filt[:,:,6] = 0.0
# specF_filt[:,:,9] = 0.0 
# specF_filt[:,:,10] = 0.0

# # Set modes to zero:          1 mode
# specF_filt[:,:,0] = 0.0
# # specF_filt[:,:,1] = 0.0
# specF_filt[:,:,2] = 0.0
# specF_filt[:,:,3] = 0.0
# # specF_filt[:,:,4] = 0.0
# specF_filt[:,:,5] = 0.0
# specF_filt[:,:,6] = 0.0
# specF_filt[:,:,7] = 0.0
# # specF_filt[:,:,8] = 0.0
# specF_filt[:,:,9] = 0.0
# specF_filt[:,:,10] = 0.0

# Set modes to zero:          1 mode in streamwise, no y forcing
specF_filt[:,:,0] = 0.0
# specF_filt[:,:,1] = 0.0
specF_filt[:,:,2] = 0.0
specF_filt[:,:,3] = 0.0
specF_filt[:,:,4] = 0.0
specF_filt[:,:,5] = 0.0
specF_filt[:,:,6] = 0.0
specF_filt[:,:,7] = 0.0
# specF_filt[:,:,8] = 0.0
specF_filt[:,:,9] = 0.0
specF_filt[:,:,10] = 0.0

# Rotate and keep only Fxp (in the streamwise)
theta = 24.0/180*np.pi
Fxp = specF_filt[:,:,1]*np.cos(theta)+specF_filt[:,:,8]*np.sin(theta)
# Project back
specF_filt[:,:,1] = np.cos(theta)*(Fxp)
specF_filt[:,:,8] = np.sin(theta)*(Fxp)

# ########################
# ### Gaussian forcing ###
# ########################
# # Base the other modes on the 1D_MWF gaussian forcing:
# # Amplitude of forces
# z = np.linspace(0,64,360)
# Fx0 = 0.035
# Fy0 = 0.002*20
# Fz0 = 0.008
# beta=np.pi/2
# def gauss(x,m,s):
#     return np.exp(-(x-m)**2/s**2/2)

# # Define Fx
# Fx = Fx0*gauss(z,47,3.5)
# # FFT
# Fxh =  np.fft.rfft(Fx)/len(Fx)
# Fxh = Fxh[:nz_c] # Truncate to desired resolution

# # Define Fy and Fz
# Fy = Fy0*np.gradient(gauss(z,46.5,4),z)/2.5
# Fz = Fz0*gauss(z,46.5,4)

# plt.plot(z,Fx)
# plt.title("Fx")
# plt.show()
# plt.plot(z,Fy)
# plt.title("Fy")
# plt.show()
# plt.plot(z,Fz)
# plt.title("Fz")
# plt.show()

# # FFT
# Fyh =  np.fft.rfft(Fy)/len(Fy)
# Fyh = Fyh[:nz_c]
# Fzh =  np.fft.rfft(Fz)/len(Fz)
# Fzh = Fzh[:nz_c]

# # Now add to modes 1,4,8
# specF_filt[:,0,1] = Fxh
# specF_filt[:,0,4] = Fyh
# specF_filt[:,0,8] = Fzh

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
