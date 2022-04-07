import numpy as np
from scipy.io import netcdf
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob as glob

def read_cdf(filename,intype='mpt'):
    """
    Loads mean-pol-tor formulation from file
    """
    with netcdf.NetCDFFile(filename,'r') as f:
        alpha = np.copy(f.alpha)
        gamma = np.copy(f.gamma)
        Re = np.copy(f.Re)
        field = np.copy(f.variables[intype].data)
        time = np.copy(f.t)
    Lx = 2*np.pi/alpha
    Lz = 2*np.pi/gamma
    return time,Re,Lx,Lz,field

def mpt2sp(mp,Lx,Lz):
    """
    Converts mean-poloidal-toroidal to u,v,w
    """
    mpRe = np.squeeze(mp[0,:,:,:])
    mpIm = np.squeeze(mp[1,:,:,:])
    
    K0 = int((mpRe.shape[2]+1)/2)
    MT = mpRe.shape[1] # 2*(MM-1), where MM is the resolution in x (e.g. Nx = 512)
    NT=  mpRe.shape[0]
    
    if MT%2==0:
        MM=int(MT/2)
    else:
        print("ERROR i_MM not divisible by two!")
        MM=int((MT+1)/2)
        
    # First derivative operator x
    ad_m1 = np.fft.fftfreq(MT,d=Lx/(MT*2*np.pi))
    ad_m2 = -ad_m1**2
    
    # Same for z
    ad_n1 = np.arange(NT)*(2*np.pi/Lz)
    ad_n2 = (np.arange(NT)*(2*np.pi/Lz))**2
    
    sRe = np.zeros((NT,MT,3*K0-1))
    sIm = np.zeros((NT,MT,3*K0-1))
    d_beta = np.pi/2
    
    klist = np.arange(K0)
    ad_k1 = (-1)**(klist) * klist * d_beta
    
    adN,adM,adK = np.meshgrid(ad_n1,ad_m1,ad_k1)
    adN2,adM2,_ = np.meshgrid(ad_n2,ad_m2,ad_k1)
    adN = np.transpose(adN,(1,0,2))
    adM = np.transpose(adM,(1,0,2))
    adK = np.transpose(adK,(1,0,2))
    adN2 = np.transpose(adN2,(1,0,2))
    adM2 = np.transpose(adM2,(1,0,2))
        
    sRe[:,:,K0-1:2*K0-1] = -mpRe[:,:,K0-1:2*K0-1]*(adM2+adN2)
    sIm[:,:,K0-1:2*K0-1] = -mpIm[:,:,K0-1:2*K0-1]*(adM2+adN2)

    sRe[:,:,:K0] = -(-mpIm[:,:,:K0]*adN +mpIm[:,:,K0-1:2*K0-1]*adK*adM)
    sIm[:,:,:K0] = (-mpRe[:,:,:K0]*adN +mpRe[:,:,K0-1:2*K0-1]*adK*adM)

    sRe[:,:,2*K0-1:3*K0-1] = -(mpIm[:,:,:K0]*adM +mpIm[:,:,K0-1:2*K0-1]*adK*adN)
    sIm[:,:,2*K0-1:3*K0-1] = (mpRe[:,:,:K0]*adM +mpRe[:,:,K0-1:2*K0-1]*adK*adN)
    
    # Constrain mean fields to be zero
    sRe[0,0,0]=sIm[0,0,0] = 0
    sRe[0,0,2*K0-1] = sIm[0,0,2*K0-1] = 0
    sRe[0,0,1:K0]= mpRe[0,0,1:K0]
    sIm[0,0,1:K0] = 0
    sRe[0,0,2*K0:3*K0-1] = mpRe[0,0,K0:2*K0-1]
    sIm[0,0,2*K0:3*K0-1] = 0
    
    return sRe+1j*sIm

def GridUy(filename,vel_field,yloc,intype='mpt'):
    """
    Produces a velocity field on a y-plane.
    
    Location (-1 to 1)
    """
    if vel_field=='U':
        ymode=0
    elif vel_field == 'W':
        ymode = 7
    elif vel_field == 'V':
        ymode = 4
    else:
        print("Not valid vel_field, returning")
        return
    
    if 'state' in filename:
        # Load file
        time,Re, Lx,Lz,mpt = read_cdf(filename,intype= intype)
        # Convert to u,v,w
        s = mpt2sp(mpt,Lx,Lz)
    else:
        # Load file
        time,Re,Lx,Lz,spec = read_cdf(filename,intype=intype)
        s = spec[0,:,:,:] + 1j*spec[1,:,:,:]
    
    # Extract velocity field at y
    if (vel_field=='U' or vel_field=='W'):
        spy = np.squeeze(
            s[:,:,ymode] 
            + np.sin(np.pi*yloc/2)*s[:,:,ymode+1] 
            + np.cos(np.pi*yloc)*s[:,:,ymode+2]
            + np.sin(3*np.pi*yloc/2)*s[:,:,ymode+3]
                        )
    else:
        spy = np.squeeze(
            + np.cos(np.pi*yloc/2)*s[:,:,ymode] 
            + np.sin(np.pi*yloc)*s[:,:,ymode+1]
            + np.cos(3*np.pi*yloc/2)*s[:,:,ymode+2]
                        )    
    
    # Transform in x
    nnx =int(np.ceil(spy.shape[1]/2))
    nnz =int(np.ceil(spy.shape[0]))
    spy_pad = np.hstack((spy[:,:nnx],np.zeros((spy.shape[0],1)),spy[:,-nnx+1:]))
    gxsz = np.fft.ifft(spy_pad,axis=1)*spy_pad.shape[1]
    
    # Transform in z
    U = np.fft.irfft(gxsz,n=2*nnz,axis=0)*2*nnz
    
    # Make grid
    Nz,Nx = U.shape
    z = np.arange(Nz)*Lz/Nz
    x = np.arange(Nx)*Lx/Nx
    X,Z = np.meshgrid(x,z)
    
    return time,Re,Lx,Lz,X,Z,U

def GridReStress_xavg(filename,intype,Ny=None):
    """
    Produces the Reynolds averages and stresses in real space (z,y).
   
    GridReStress_xavg(filename, intype,Ny=None).
    If filename is even, intypes are: umean,wmean,uumean,uwmean,wwmean,vvmean

    If filename is odd, intypes are: vmean,wvmean,uvmean
    
    By default, the vertical resolution is based on K0, but if you specify Ny=16, it will make 16 y points in the vertical.

    Returns: time,Re,Lx,Lz,Y,Z,r 
    """
    # Load file
    time,Re,Lx,Lz,field = read_cdf(filename,intype=intype)
    
    NN = field.shape[1]
    K  = field.shape[2]
    
    s = field[0] + 1j*field[1]
    
    Ly = 2

    if Ny==None:
        Ny = (K-1)*2
        y = Ly*np.arange(Ny)/(Ny-1) - Ly/2
    else:
        y = np.linspace(-1.,1,Ny)
        
    ry = np.zeros((NN,Ny),dtype=complex)
    # convert to real space in y.
    if (K%2)==0:
        for i in range(len(y)):
            ry[:,i] = s[:,0] + np.sin(np.pi*y[i]/2)*s[:,1] + np.cos(np.pi*y[i])*s[:,2] + np.sin(3*np.pi*y[i]/2)*s[:,3]
    else:
        for i in range(len(y)):
            ry[:,i] = np.cos(np.pi*y[i]/2)*s[:,0] + np.sin(np.pi*y[i])*s[:,1] + np.cos(3*np.pi*y[i]/2)*s[:,2]

    # Transform in z
    nnz =int(np.ceil(ry.shape[0]))

    # Transform in z
    r = np.fft.irfft(ry,n=2*nnz,axis=0)*2*nnz

    # Make grid
    Nz,_ = r.shape
    z = np.arange(Nz)*Lz/Nz
    Z,Y = np.meshgrid(z,y,indexing='ij')

    
    return time,Re,Lx,Lz,Y,Z,r


def GridReForces_xavg(run,filenum,direction,Ny=None):
    """
    Produces the forces due to (x-avged) Reynolds stresses in real space (z,y).
   
    GridReForces_xavg(run,filenum, direction,Ny=None). Direction = 'x', 'y', or 'z'
    
    By default, the vertical resolution is based on K0, but if you specify Ny=16, it will make 16 y points in the vertical.

    Returns: time,Re,Lx,Lz,Y,Z,F 
    """
    filenum = str(filenum).zfill(4)
    restress_even = sorted(glob.glob(run+'/restress_even'+filenum+'.cdf.dat'))
    restress_odd = sorted(glob.glob(run+'/restress_odd'+filenum+'.cdf.dat'))
    
    if direction == 'x':
        parity = {'uvmean':'odd','uwmean':'even'}
    elif direction == 'y':
        parity = {'vvmean':'even','wvmean':'odd'}
    elif direction == 'z':
        parity = {'wvmean':'odd','wwmean':'even'}
    else:
        print("Please specify a direction, x, y, or z.")
        return
    
    for ii,field in enumerate(parity):
        # Filename:
        if parity[field]=='odd':
            filename = restress_odd[0]
        elif parity[field]=='even':
            filename = restress_even[0]

        # Load file:
        time,Re,Lx,Lz,dat = read_cdf(filename,intype=field)

        NN = dat.shape[1]
        K  = dat.shape[2]

        s = dat[0] + 1j*dat[1]

        # Prepare z derivatives:
        ad_n1 = np.arange(NN)*(2*np.pi/Lz)
        adN,adK = np.meshgrid(ad_n1,s.shape[1], indexing='ij')

        if ii==1: # Take z-derivative
            s *= adN*1j

        # And finally transform to real space:
        Ly = 2

        if Ny==None:
            Ny = (K-1)*2
            y = Ly*np.arange(Ny)/(Ny-1) - Ly/2
        else:
            y = np.linspace(-1.,1,Ny)

        ry = np.zeros((NN,Ny),dtype=complex)
        if ii==0:
            # convert to real space in y, but taking the derivative in y
            if (K%2)==0:
                for i in range(len(y)):
                    ry[:,i] = np.pi/2*np.cos(np.pi*y[i]/2)*s[:,1] - np.pi*np.sin(np.pi*y[i])*s[:,2] + 3*np.pi/2*np.cos(3*np.pi*y[i]/2)*s[:,3]
            else:
                for i in range(len(y)):
                    ry[:,i] = -np.pi/2*np.sin(np.pi*y[i]/2)*s[:,0] + np.pi*np.cos(np.pi*y[i])*s[:,1] - 3*np.pi/2*np.sin(3*np.pi*y[i]/2)*s[:,2]
        else:    
            # convert to real space in y.
            if (K%2)==0:
                for i in range(len(y)):
                    ry[:,i] = s[:,0] + np.sin(np.pi*y[i]/2)*s[:,1] + np.cos(np.pi*y[i])*s[:,2] + np.sin(3*np.pi*y[i]/2)*s[:,3]
            else:
                for i in range(len(y)):
                    ry[:,i] = np.cos(np.pi*y[i]/2)*s[:,0] + np.sin(np.pi*y[i])*s[:,1] + np.cos(3*np.pi*y[i]/2)*s[:,2]

        # Transform in z
        nnz =int(np.ceil(ry.shape[0]))

        # Transform in z
        r = np.fft.irfft(ry,n=2*nnz,axis=0)*2*nnz

        # Make grid
        Nz,_ = r.shape
        z = np.arange(Nz)*Lz/Nz
        Z,Y = np.meshgrid(z,y,indexing='ij')

        if ii==0:
            final = np.copy(r)
        if ii==1:
            final += np.copy(r)
    
    return time,Re,Lx,Lz,Y,Z,final

def SpecReForces_xavg(run,filenum,direction):
    """
    Produces the forces due to Reynolds stresses in real space (z) but keeps the vertical direction spectral.
   
    SpecReForces_xavg(run,filenum, direction). Direction = 'x', 'y', or 'z'
    
    Returns: time,Re,Lx,Lz,Y,Z,F 
    """
    filenum = str(filenum).zfill(4)
    restress_even = sorted(glob.glob(run+'/restress_even'+filenum+'.cdf.dat'))
    restress_odd = sorted(glob.glob(run+'/restress_odd'+filenum+'.cdf.dat'))
    
    if direction == 'x':
        parity = {'uvmean':'odd','uwmean':'even'}
    elif direction == 'y':
        parity = {'vvmean':'even','wvmean':'odd'}
    elif direction == 'z':
        parity = {'wvmean':'odd','wwmean':'even'}
    else:
        print("Please specify a direction, x, y, or z.")
        return
    
    for ii,field in enumerate(parity):
        # Filename:
        if parity[field]=='odd':
            one = -1
            filename = restress_odd[0]
        elif parity[field]=='even':
            one = 1
            filename = restress_even[0]

        # Load file:
        time,Re,Lx,Lz,dat = read_cdf(filename,intype=field)

        NN = dat.shape[1]
        K  = dat.shape[2]

        s = dat[0] + 1j*dat[1]
        
        if (K%2)!=0: # Add the zero mode to y.
            s=np.hstack([np.zeros((NN,1),dtype=complex),s])
            K += 1

        # Prepare z derivatives:
        ad_n1 = np.arange(NN)*(2*np.pi/Lz)
        klist = np.arange(K)
        ad_k1 = (-1)**(klist+1)*klist*(np.pi/2)
        adN,adK = np.meshgrid(ad_n1,ad_k1, indexing='ij')

        if ii==1: # Take z-derivative
            s *= adN*1j
        elif ii==0: # Take y-derivative
            s *= one*adK

        # Transform in z
        nnz =int(np.ceil(s.shape[0]))

        # Transform in z
        r = np.fft.irfft(s,n=2*nnz,axis=0)*2*nnz
            
        # Make grid
        Nz,_ = r.shape
        z = np.arange(Nz)*Lz/Nz
#         Z,Y = np.meshgrid(z,y,indexing='ij')

        if ii==0:
            final = np.copy(r)
        if ii==1:
            final += np.copy(r)
    
    return time,Re,Lx,Lz,z,final

def GridXavg(filename,vel_field,Ny=None,intype='mpt'):
    """
    Produces an x-averaged velocity field on a y-z plane. 
    
    GridXavg(filename,vel_field,Ny=None). Possible vel_fields are U,V,W
    """
    if vel_field=='U':
        ymode=0
        K = 4
    elif vel_field == 'W':
        ymode = 7
        K = 4
    elif vel_field == 'V':
        ymode = 4
        K=3
    else:
        print("Not valid vel_field, returning")
        return
    
    if 'state' in filename:
        # Load file
        time,Re, Lx,Lz,mpt = read_cdf(filename,intype= intype)
        # Convert to u,v,w
        s = mpt2sp(mpt,Lx,Lz)
    else:
        # Load file
        time,Re,Lx,Lz,spec = read_cdf(filename,intype=intype)
        s = spec[0,:,:,:] + 1j*spec[1,:,:,:]
    
    Ly = 2
    
    if Ny==None:
        Ny = (K-1)*2
        y = Ly*np.arange(Ny)/(Ny-1) - Ly/2
    else:
        y = np.linspace(-1.,1,Ny)
    
    ry = np.zeros((s.shape[0],Ny),dtype=complex)
    # Convert the horizontal mean to real space in y 
    if (vel_field=='U' or vel_field=='W'):
        for ii in range(len(y)):
            ry[:,ii] = s[:,0,ymode] + np.sin(np.pi*y[ii]/2)*s[:,0,ymode+1] + np.cos(np.pi*y[ii])*s[:,0,ymode+2]+ np.sin(3*np.pi*y[ii]/2)*s[:,0,ymode+3]
    else:
        for ii in range(len(y)):
            ry[:,ii] = np.cos(np.pi*y[ii]/2)*s[:,0,ymode] + np.sin(np.pi*y[ii])*s[:,0,ymode+1] + np.cos(3*np.pi*y[ii]/2)*s[:,0,ymode+2]
    
    # Transform in z
    nnz =int(np.ceil(ry.shape[0]))

    # Transform in z
    r = np.fft.irfft(ry,n=2*nnz,axis=0)*2*nnz

    # Make grid
    Nz,_ = r.shape
    z = np.arange(Nz)*Lz/Nz
    Z,Y = np.meshgrid(z,y,indexing='ij')
    
    return time,Re,Lx,Lz,Y,Z,r

def plot_reforces(run,filenum,direction,figwidth = 12,aspect_factor=3,Ny=None):
    """
    For help choosing arguments, look at help(GridReForces_xavg).
    """
    time,Re,Lx,Lz,Y,Z,r = GridReForces_xavg(run,filenum,direction,Ny=Ny)
    
    rbar = np.max([-np.min(r),np.max(r)])
    
    plt.figure(figsize=(figwidth,figwidth/3))
    plt.title("Re = %.2e, Lx = %.2f, Lz = %.2f, \n Direction: %s, time= %.4f" % (Re,Lx,Lz,direction,time))
    if Ny==None:
        plt.imshow(r.T,vmin=-rbar,vmax=rbar,cmap='bwr',aspect=(np.shape(Z)[0]/np.shape(Z)[1]/aspect_factor))
    else:
        plt.pcolor(Z,Y,r,vmin=-rbar,vmax=rbar,cmap='bwr')
#         plt.gca().set_aspect(np.shape(Z)[0]/np.shape(Z)[1]/aspect_factor)
    plt.colorbar()
    plt.ylabel(r'$y$')
    plt.xlabel(r'$z$')
    plt.show()

def plot_restress(filename,field,figwidth = 12,aspect_factor=3,Ny=None):
    """
    For help choosing filenames and fields, look at help(GridReStress_xavg).
    """
    time,Re,Lx,Lz,Y,Z,r = GridReStress_xavg(filename,field,Ny=Ny)
    
    rbar = np.max([-np.min(r),np.max(r)])
    
    plt.figure(figsize=(figwidth,figwidth/3))
    plt.title("Re = %.2e, Lx = %.2f, Lz = %.2f, \n Field: %s, time= %.4f" % (Re,Lx,Lz,field,time))
    if Ny==None:
        plt.imshow(r.T,vmin=-rbar,vmax=rbar,cmap='bwr',aspect=(np.shape(Z)[0]/np.shape(Z)[1]/aspect_factor))
    else:
        plt.pcolor(Z,Y,r,vmin=-rbar,vmax=rbar,cmap='bwr')
#         plt.gca().set_aspect(np.shape(Z)[0]/np.shape(Z)[1]/aspect_factor)
    plt.colorbar()
    plt.ylabel(r'$y$')
    plt.xlabel(r'$z$')
    plt.show()
    
def plot_XZ(filename,field,yloc,intype='mpt',figwidth=8):
    """
    If filename is state****.cdf.dat:
        Plots field at yloc from filename.
        
    If filename is restress_2d:
        Plots the 2D Reynolds average or Reynolds stresses at a given y location ('yloc'). 
        Must specify an 'intype'. Options include 'umean_2d', which includes (umean_2d,vmean_2d,wmean_2d), 'remean1_2d', which includes (uumean,uvmean,uwmean), and 'remean2_2d', which includes (vvmean,wvmean,wwmean). To choose between the three options, specify 'U', 'V', or 'W' for 'field', for positions 1, 2, and 3 in the lists above. 
    """
    time,Re,Lx,Lz,X,Z,U = GridUy(filename,field,yloc,intype=intype)
    
    cscale=np.max((np.abs(np.min(U)),np.max(U)))
    
    plt.figure(figsize=(figwidth,figwidth*(Lz/Lx)))
    plt.title("Re = %.2e, Lx = %.2f, Lz = %.2f, \n Field: %s, yloc = %s, time= %.4f" % (Re,Lx,Lz,field,yloc,time))
    plt.pcolormesh(X,Z,U,vmin=-cscale,vmax=cscale,cmap=cm.bwr)
    plt.colorbar()
    ax=plt.gca()
    ax.set_aspect(1)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$z$')
    plt.show()
    
def plot_YZ_xavg(filename,field,intype='mpt',figwidth = 12,aspect_factor=3,Ny=None):
    """
    If filename is state****.cdf.dat:
        Plots x-averaged field from filename, with a given Ny resolution (default = None, which means the choice is based on the spectral resolution).
        
    If filename is restress_2d:
        Plots the x-averaged Reynolds average or Reynolds stresses. 
        Must specify an 'intype'. Options include 'umean_2d', which includes (umean_2d,vmean_2d,wmean_2d), 'remean1_2d', which includes (uumean,uvmean,uwmean), and 'remean2_2d', which includes (vvmean,wvmean,wwmean). To choose between the three options, specify 'U', 'V', or 'W' for 'field', for positions 1, 2, and 3 in the lists above. 
    
    For help choosing filenames and fields, look at help(GridXavg).
    """
    time,Re,Lx,Lz,Y,Z,r = GridXavg(filename,field,Ny=Ny,intype=intype)
    
    rbar = np.max([-np.min(r),np.max(r)])
    
    plt.figure(figsize=(figwidth,figwidth/3))
    plt.title("Re = %.2e, Lx = %.2f, Lz = %.2f, \n Field: %s, time= %.4f" % (Re,Lx,Lz,field,time))
    if Ny==None:
        plt.imshow(r.T,vmin=-rbar,vmax=rbar,cmap='bwr',aspect=(np.shape(Z)[0]/np.shape(Z)[1]/aspect_factor))
    else:
        plt.pcolor(Z,Y,r,vmin=-rbar,vmax=rbar,cmap='bwr')
#         plt.gca().set_aspect(np.shape(Z)[0]/np.shape(Z)[1]/aspect_factor)
    plt.colorbar()
    plt.ylabel(r'$y$')
    plt.xlabel(r'$z$')
    plt.show()
    
def plot_turb(filename,figwidth=8):
    time,Re,Lx,Lz,turb = read_cdf(filename,intype='uvw')
    turb_tot = np.sum(turb,axis=2)
    x = np.linspace(0,Lx,turb_tot.shape[0])
    z = np.linspace(0,Lz,turb_tot.shape[1])
    Z,X = np.meshgrid(z,x)
    
    cscale=np.max(turb_tot)
    
    plt.figure(figsize=(figwidth,figwidth*(Lz/Lx)))
    plt.title("TKE, Re = %.2e, Lx = %.2f, Lz = %.2f, time= %.4f" % (Re,Lx,Lz,time),fontsize=15)
    plt.pcolormesh(X,Z,turb_tot,vmin=0,vmax=cscale,cmap=cm.Reds)
    plt.colorbar()
    ax=plt.gca()
    ax.set_aspect(1)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$z$')
    plt.show()
