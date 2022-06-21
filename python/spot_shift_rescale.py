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
odir = '../'
# Rescaling factor
rfactor = 4
# Amplitfying factor
amp = 10
print("scaling down by a factor of %s" % rfactor)
# Spectral truncation mode numbers
MM_new = int(2048/rfactor)
M_new = int(2*(MM_new-1))
N_new = int(2048/rfactor) 

# Name for  saved file
name = 'spot_IC_{}_{}_amp_{}.bak'.format(MM_new,N_new,amp)
# Name of original file to use
spot_IC = idir+'spot_IC.bak'

################################################
### Finding max, center of mass and shifting ###
################################################
print("Performing first shift...")
# First approximation for centering.
# Load data
yloc=0.0
vel_field='U'
time,Re,Lx,Lz,s = aa.load_spec(spot_IC,intype='mpt')
time,Re,Lx,Lz,X,Z,U = aa.GridUy(time,Re,Lx,Lz,s,vel_field,yloc)
del s
# Location of shift (indices)
inds = np.unravel_index(np.abs(U).argmax(), np.abs(U).shape)
xmax = X[inds]
zmax = Z[inds]
del X,Z,U

# Center the spot
xfinal = Lx/2
zfinal = Lz/2

xshift = xmax - xfinal
zshift = zmax - zfinal

print("(xfinal,zfinal) = ({},{})".format(xfinal,zfinal))

print("(xshift,zshift) = ({},{})".format(xshift,zshift))

# Apply shift to spot_IC
# Load file
time,Re, Lx,Lz,mp = aa.read_cdf(spot_IC,intype= 'mpt')

# Write in complex form
mpt = np.squeeze(mp[0,:,:,:]) + 1j*np.squeeze(mp[1,:,:,:])

K0 = int((mpt.shape[2]+1)/2)
MT = mpt.shape[1] # 2*(MM-1) = i_M, where MM is the resolution in x (e.g. Nx = 512). Goes from 0 to i_M1 = i_M -1 = 2*(i_MM-1)-1
NT=  mpt.shape[0]

if MT%2==0:
    MM=int(MT/2)+1
else:
    print("ERROR i_MM not divisible by two!")
    MM=np.nan

# Prepare x-wavenumbers
ad_m1 = 1j*np.fft.fftfreq(MT,d=Lx/(MT*2*np.pi))

# Prepare z-wavenumbers
ad_n1 = 1j*np.arange(NT)*(2*np.pi/Lz)

adN,adM,adK = np.meshgrid(ad_n1,ad_m1,np.arange(mpt.shape[2]),indexing='ij')

# Make shift matrix, exp(i*kx*x + i*kz*z)
mat = np.exp(adM*xshift + adN*zshift)

# Apply shift:
mpt *= mat

del mp,mat

# Now that we have an approximately shifted spot, let's fine-tune the shift by calculating the Center of Mass and making the shift more accurate
print("Performing second shift...")
# Transform in x
KK = mpt.shape[2] 
nnx =int(np.ceil(mpt.shape[1]/2)) # MM
nnz =int(np.ceil(mpt.shape[0]))
spy_pad = np.hstack((mpt[:,:nnx,:],np.zeros((mpt.shape[0],1,KK)),mpt[:,-nnx+1:,:]))
gxsz = np.fft.ifft(spy_pad,axis=1)*spy_pad.shape[1]

# Transform in z
U = np.fft.irfft(gxsz,n=2*nnz,axis=0)*2*nnz

Nz,Nx = U[:,:,0].shape
z = np.arange(Nz)*Lz/Nz
x = np.arange(Nx)*Lx/Nx
X,Z = np.meshgrid(x,z)

del gxsz,spy_pad

xcom = np.sum(np.abs(U[:,:,0])*X)/np.sum(np.abs(U[:,:,0]))
zcom = np.sum(np.abs(U[:,:,0])*Z)/np.sum(np.abs(U[:,:,0]))

# Center the spot
xfinal = Lx/2
zfinal = Lz/2

xshift = xcom - xfinal
zshift = zcom - zfinal

print("Second shift: (xfinal,zfinal) = ({},{})".format(xfinal,zfinal))

print("Second shift: (xshift,zshift) = ({},{})".format(xshift,zshift))

# Apply shift to spot_IC
# Make shift matrix, exp(i*kx*x + i*kz*z)
mat = np.exp(adM*xshift + adN*zshift)

# Apply shift:
mpt *= mat

del mat,adM,adN,adK

##################
### Truncation ###
##################
print("Truncating in real space...")
# Okay, we have shifted the spot. Next we FFT
# Transform in x
KK = mpt.shape[2] 
nnx =int(np.ceil(mpt.shape[1]/2)) # MM
nnz =int(np.ceil(mpt.shape[0]))
spy_pad = np.hstack((mpt[:,:nnx,:],np.zeros((mpt.shape[0],1,KK)),mpt[:,-nnx+1:,:]))
gxsz = np.fft.ifft(spy_pad,axis=1)*spy_pad.shape[1]

# Transform in z
U = np.fft.irfft(gxsz,n=2*nnz,axis=0)*2*nnz

del gxsz

# Truncate in real space
Nz,Nx,Ny = U.shape
U_trun = U[int(Nz*(rfactor - 1)/(2*rfactor)):-int(Nz*(rfactor - 1)/(2*rfactor)),int(Nz*(rfactor - 1)/(2*rfactor)):-int(Nz*(rfactor - 1)/(2*rfactor)),:]

# # Now we need to apply some kind of windowing function
# z = np.arange(U_trun.shape[0])
# x = np.arange(U_trun.shape[1])
# Z,X,Y = np.meshgrid(z,x,np.zeros(7),indexing='ij')

# U_trun *= np.exp(-((X-Lx/rfactor)**2 + (Z-Lz/rfactor)**2)/(2048/2/rfactor/4)**2)

Wz = np.bartlett(U_trun.shape[0])
Wx = np.bartlett(U_trun.shape[1])
WZ,WX,WY = np.meshgrid(Wz,Wx,np.ones(7),indexing='ij')

U_trun *= WZ*WX*amp
# There will still be some error, but hopefully this minimizes it... It makes the 'jump' about two orders of magnitude smaller than the magnitude at the spot. We'll see if it's good enough

# Transform back
# Transform in z.
gxsz = np.fft.rfft(U_trun,axis=0)
# Truncate extra mode
gxsz = gxsz[:-1,:,:]
gxsz *= 1/(2*nnz)
# Transform in x
mpt_2 = np.fft.fft(gxsz,axis=1)
mpt_2 *= 1/U.shape[1]
# Set zero modes to zero!
mpt_2[0,0,0] = 0.0

# Finally, reduce the resolution
mpt_final = np.zeros((N_new,M_new,KK),dtype=complex)
mpt_final[:,:MM_new,:] = mpt_2[:N_new,:MM_new,:]
mpt_final[:,MM_new:,:] = mpt_2[:N_new,-MM_new+2:,:]

del mpt_2,U,U_trun,gxsz

# Define new domain sizes for feeding into the cdf file
Lx *= 1/rfactor
Lz *= 1/rfactor

##############
### Saving ###
##############
print("Saving...")
with netcdf.NetCDFFile(odir+name,'w') as f:
    # Create dimensions
    f.createDimension('N', N_new) #NN
    f.createDimension('M', M_new) #M
    f.createDimension('K',KK) #K
    f.createDimension('ReIm',2)
    
    # define attributes
    f.t = 0.0
    f.t_step = 0
    f.Re = 80.0
    f.alpha = 2*np.pi/Lx
    f.gamma = 2*np.pi/Lz
    
    # Define variable
    mpt_in = f.createVariable('mpt', 'd', ('ReIm','N','M','K'))
    mpt_in.K=KK
    mpt_in.MM=MM_new
    mpt_in.NN=N_new
    mpt_in[0,:,:,:] = np.real(mpt_final)
    mpt_in[1,:,:,:] = np.imag(mpt_final)
    
print("Done! Saved %s in %s" % (name,odir))
