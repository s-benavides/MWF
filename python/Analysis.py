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
        try:
            field = np.copy(f.variables[intype].data)
        except:
            print("Not a valid 'intype'! These are the possible options:")
            print(f.variables.keys())
        time = np.copy(f.t)
    Lx = 2*np.pi/alpha
    Lz = 2*np.pi/gamma
    return time,Re,Lx,Lz,field

def mpt2sp(mp,Lx,Lz):
    """
    Converts mean-poloidal-toroidal to u,v,w
    """
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
    
    # Prepare x-derivative
    ad_m1 = 1j*np.fft.fftfreq(MT,d=Lx/(MT*2*np.pi))
    ad_m2 = ad_m1**2
    
    # Prepare z-derivative
    ad_n1 = 1j*np.arange(NT)*(2*np.pi/Lz)
    ad_n2 = ad_n1**2
    # Prepare y derivatives:
    klist = np.arange(K0)
    ad_k1 = (-1)**(klist)*klist*(np.pi/2) # Include ky = 0 mode.

    adN,adM,adK = np.meshgrid(ad_n1,ad_m1,ad_k1,indexing='ij')
    adN2,adM2,_ = np.meshgrid(ad_n2,ad_m2,ad_k1,indexing='ij')
    
    ## Spec array
    s = np.zeros((NT,MT,3*K0-1),dtype=complex)
    
    # u
    s[:,:,:K0] = -adN*mpt[:,:,:K0] + adM*adK*mpt[:,:,K0-1:]
    # v
    s[:,:,K0:2*K0-1] = -(adM2[:,:,1:]+adN2[:,:,1:])*mpt[:,:,K0:]
    # w
    s[:,:,2*K0-1:] = adM*mpt[:,:,:K0] + adK*adN*mpt[:,:,K0-1:] 
    
    # Constrain mean fields to be zero
    # u
    s[0,0,0] = 0.0
    # w 
    s[0,0,2*K0-1] = 0.0
    
    # Add f(y) and g(y)
    s[0,0,1:K0]= np.real(mpt[0,0,1:K0])
    s[0,0,2*K0:3*K0-1] = np.real(mpt[0,0,K0:])
    
    return s

def load_spec(filename,intype='mpt',Retype=None):    
    """
    Load a field into spec type. Either read a state file, restress file, or calculate the spectral forces if 'intype'==SpecReForces.
    
    Must specify an 'intype'. Options include 'umean_2d', which includes (umean_2d,vmean_2d,wmean_2d), 'remean1_2d', which includes (uumean,uvmean,uwmean), and 'remean2_2d', which includes (vvmean,wvmean,wwmean). If plotting forces, then intype = 'SpecReForces', you must specify Retype='2d' or 'filt', and field ('U','V', or 'W') determine the direction of the force. To choose between the three options, specify 'U', 'V', or 'W' for 'field', for positions 1, 2, and 3 in the lists above. 
    """
    if 'state' in filename:
        # Load file
        time,Re, Lx,Lz,mpt = read_cdf(filename,intype= intype)
        # Convert to u,v,w
        s = mpt2sp(mpt,Lx,Lz)
    elif intype=='SpecReForces':
        time,Re,Lx,Lz,s = SpecReForces(filename,Retype=Retype)
    else:
        # Load file
        time,Re,Lx,Lz,spec = read_cdf(filename,intype=intype)
        s = spec[0,:,:,:] + 1j*spec[1,:,:,:]
        
    return time,Re,Lx,Lz,s

def GridUy(time,Re,Lx,Lz,s,vel_field,yloc):
    """
    Produces a velocity field on a y-plane, given a spec field (use 'load_spec' function to load and generate the field from saved files, or create your own).
    
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

def SpecReForces_xavg(run,filenum,grid=False):
    """
    Produces the forces due to Reynolds stresses (calculated from the x-averaged stresses 'restress_even' or 'restress_odd') in spectral space if 'grid=False' and in real space (z) if 'grid=True' but keeps the vertical direction spectral.
   
    SpecReForces_xavg(run,filenum).
    
    Returns: time,Re,Lx,Lz,z,out 
    """
    filenum = str(filenum).zfill(4)
    restress_even = sorted(glob.glob(run+'/restress_even'+filenum+'.cdf.dat'))
    restress_odd = sorted(glob.glob(run+'/restress_odd'+filenum+'.cdf.dat'))
    
    for direction in ['x','y','z']:
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
            if parity[field]=='even':
                KK=int(3*K-1)
                
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
            
            if grid:
                # Transform in z
                nnz =int(np.ceil(s.shape[0]))

                # Transform in z
                s = np.fft.irfft(s,n=2*nnz,axis=0)*2*nnz

                # Make grid
                Nz,_ = s.shape
                z = np.arange(Nz)*Lz/Nz
            else:
                z = np.nan

            if ii==0:
                final = np.copy(s)
            if ii==1:
                final += np.copy(s)
        
        if direction=='x':
            out = np.zeros((s.shape[0],KK),dtype=complex)
            out[:,0:4] = np.copy(final)
        if direction=='y':
            # If odd mode, remove the 'zero' mode used in calculations.
            out[:,4:7] = np.copy(final[:,1:])
        elif direction=='z':
            out[:,7:11] = np.copy(final)
    
    return time,Re,Lx,Lz,z,out

def SpecReForces(filename,Retype='2d'):
    """
    Produces the forces due to Reynolds stresses in spectral space space.
   
    SpecReForces(filename,Retype='2d'). Outputs the forces in spectral space, consistent with 'spec' type in MWF. In other words, 11 vertical modes, Fx in the first four, Fy in the next three, and Fz in the final four.
    
    Retype could also be 'filt'.
    
    Returns: time,Re,Lx,Lz,F 
    """
    for direction in ['x','y','z']:
        if direction == 'x':
            recomp = {'uu':['even','remean1','U'],'uv':['odd','remean1','V'],'uw':['even','remean1','W']}
        elif direction == 'y':
            recomp = {'uv':['odd','remean1','V'],'vv':['even','remean2','U'],'wv':['odd','remean2','V']}
        elif direction == 'z':
            recomp = {'uw':['even','remean1','W'],'wv':['odd','remean2','V'],'ww':['even','remean2','W']}
        else:
            print("Please specify a direction: x, y, or z.")
            return

        for ii,field in enumerate(recomp):

            with netcdf.NetCDFFile(filename,'r') as f:
                KK = f.dimensions['KK']
                M = f.dimensions['M']
                N = f.dimensions['N']

            K0 = int((KK+1)/3)

            if recomp[field][2] =='U':
                ymode=0
                ypl = 4
            elif recomp[field][2] == 'W':
                ymode = 7
                ypl = 4
            elif recomp[field][2] == 'V':
                ymode = 4
                ypl = 3

            # Load file
            time,Re,Lx,Lz,spec = read_cdf(filename,intype=recomp[field][1]+'_'+Retype)
            # Isolate specific field we want.
            s = spec[0,:,:,ymode:ymode+ypl] + 1j*spec[1,:,:,ymode:ymode+ypl]

            if ypl==3: # Add the zero mode to y.
                s=np.dstack([np.zeros((N,M,1),dtype=complex),s])

            # Prepare z derivatives:
            ad_n1 = np.arange(N)*(2*np.pi/Lz)
            # Prepare x derivatives:
            ad_m1 = np.fft.fftfreq(M,d=Lx/(M*2*np.pi))
            # Prepare y derivatives:
            if recomp[field][0]=='odd':
                one = -1
            elif recomp[field][0]=='even':
                one = 1
            klist = np.arange(K0)
            ad_k1 = (-1)**(klist+1)*klist*(np.pi/2)
            adN,adM,adK = np.meshgrid(ad_n1,ad_m1,ad_k1, indexing='ij')

            if ii==0: # Take x-derivative
                s *= adM*1j
            elif ii==1: # Take y-derivative
                s *= one*adK
            elif ii==2: # Take z-derivative
                s *= adN*1j

            if ii==0:
                final = np.copy(s)
            else:
                final += np.copy(s)

        if direction=='x':
            spec_out = np.zeros((N,M,KK),dtype=complex)
            spec_out[:,:,0:4] = np.copy(final)
        if direction=='y':
            # If odd mode, remove the 'zero' mode used in calculations.
            spec_out[:,:,4:7] = np.copy(final[:,:,1:])
        elif direction=='z':
            spec_out[:,:,7:11] = np.copy(final)

    return time,Re,Lx,Lz,spec_out

def GridXavg(time,Re,Lx,Lz,s,vel_field,Ny=None):
    """
    Produces an x-averaged velocity field on a y-z plane. 
    
    GridXavg(time,Re,Lx,Lz,s,vel_field,Ny=None). Possible vel_fields are U,V,W
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
    
    dirs = {'U':'x','V':'y','W':'z'}
    
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
    
def plot_XZ(filename,field,yloc,intype='mpt',figwidth=8,Retype=None):
    """
    If filename is state****.cdf.dat:
        Plots field at yloc from filename.
        
    If filename is restress_2d:
        Plots the 2D Reynolds average or Reynolds stresses at a given y location ('yloc'). 
        Must specify an 'intype'. Options include 'umean_2d', which includes (umean_2d,vmean_2d,wmean_2d), 'remean1_2d', which includes (uumean,uvmean,uwmean), and 'remean2_2d', which includes (vvmean,wvmean,wwmean). To choose between the three options, specify 'U', 'V', or 'W' for 'field', for positions 1, 2, and 3 in the lists above. If plotting forces, then intype = 'SpecReForces', you must specify Retype='2d' or 'filt', and field ('U','V', or 'W') determine the direction of the force. To choose between the three options, specify 'U', 'V', or 'W' for 'field', for positions 1, 2, and 3 in the lists above. 
    """
    time,Re,Lx,Lz,s = load_spec(filename,intype=intype,Retype=Retype)
    time,Re,Lx,Lz,X,Z,U = GridUy(time,Re,Lx,Lz,s,field,yloc)
    
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
    
def plot_YZ_xavg(filename,field,intype='mpt',figwidth = 12,aspect_factor=3,Ny=None,Retype=None):
    """
    If filename is state****.cdf.dat:
        Plots x-averaged field from filename, with a given Ny resolution (default = None, which means the choice is based on the spectral resolution).
        
    If filename is restress_2d:
        Plots the x-averaged Reynolds average or Reynolds stresses. 
    
    For help choosing filenames and fields, look at help(GridXavg).
    """
    time,Re,Lx,Lz,s = load_spec(filename,intype=intype,Retype=Retype)
    time,Re,Lx,Lz,Y,Z,r = GridXavg(time,Re,Lx,Lz,s,field,Ny=Ny)
    
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

def spectrum_calc(spec_field,Lx,Lz):
    """
    Takes a field (of type 'spec') and outputs the power spectrum in the z direction (averaged in x) and x direction (averaged in z). Keeps the y-modes separate, so that the final spectra are two-dimensional arrays of size, e.g. (N,KK) for the z spectrum.
    
    returns kz,specz,|kx|,specx
    """
    N,M,KK = field.shape
    ad_n1 = np.arange(N)*(2*np.pi/Lz)
    ad_m1 = np.fft.fftfreq(M,d=Lx/(M*2*np.pi))
    specz = np.mean(np.abs(field),axis=1)
    
    specx = np.mean(np.abs(field),axis=0)
    kmags = np.unique(np.abs(ad_m1))
    specx_f = np.zeros(kmags.shape)
    for m in range(M):
        kmag = np.abs(ad_m1[m])
        specx_f[np.where(kmags==kmag)] += specx[m]
    
    return ad_n1,specz,kmags,specx_f
    
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
