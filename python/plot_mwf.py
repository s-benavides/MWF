import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import matplotlib.cm as cm
import Analysis as aa
from scipy.io import netcdf

plt.rcParams.update({'font.size': 20})

# Choose input directory
idir = '../'

# Choose output directory
odir = './Figures/'

# Searches through all directories in 'Data' folder (which are named after experiments) and imports the data:
#dirs = sorted(glob.glob(idir+'MFU*'))
#dirs = sorted(glob.glob(idir+'spot_test*'))
dirs = sorted(glob.glob(idir+'Lx*'))

runs = []
for file in dirs:
    run = file.split('/')[1]
    runs.append(run)

# nremove = input('How many runs to remove?')
# nremove = int(nremove)
# if nremove>0:
#     for i in range(nremove):
#         rem = input("What run to remove?")
#         runs.remove(rem)

print(runs)

nruns = input('How many runs to plot? ')
nruns = int(nruns)

runs = []
for i in range(nruns):
    run = input("Run name: ")
    runs.append(run)


###################
### Time series ###
###################
plt.figure(figsize=(10,4))
for run in runs[:]:
    t=np.genfromtxt(idir+run+'/vel_energy.dat')[:,0]
    ke=np.genfromtxt(idir+run+'/vel_energy.dat')[:,1]
    print('initial KE = %s' % ke[0])
#     print("length = %s" % len(ke))
    plt.semilogy(ke,label=run)
    
# plt.title(run)
# print('max t = %s' % np.max(t))
plt.legend(loc=(1.01,0.1),fontsize=15)
plt.tight_layout()
plt.show()

##########################
### Flow visualization ###
##########################

# figwidths=[10,10,10]
for ii,run in enumerate(runs):
    print(run)
    state_outs = sorted(glob.glob(idir+run+'/state*'))
##     print(state_outs[:])
    for state_out in np.array(state_outs)[-1:]:
#        print(state_out)
        aa.plot_XZ(state_out,'U',0.0,figwidth=10)
        aa.plot_XZ(state_out,'V',0.0,figwidth=10)
        aa.plot_XZ(state_out,'W',0.0,figwidth=10)
        aa.plot_YZ_xavg(state_out,'U',Ny=100,figwidth=24)
        aa.plot_YZ_xavg(state_out,'V',Ny=100,figwidth=24)
        aa.plot_YZ_xavg(state_out,'W',Ny=100,figwidth=24)


#     spec_in = sorted(glob.glob(idir+run+'/spec.cdf*'))
#     if len(spec_in)>0:
#         aa.plot_YZ_xavg(spec_in[-1],'U',intype='spec',figwidth=24,Ny=100)
#         aa.plot_YZ_xavg(spec_in[-1],'V',intype='spec',figwidth=24,Ny=100)
#         aa.plot_YZ_xavg(spec_in[-1],'W',intype='spec',figwidth=24,Ny=100)
    re_in = sorted(glob.glob(idir+run+'/restress_2d*'))
    if len(re_in)>0:
        aa.plot_YZ_xavg(re_in[-1],'U',intype='umean_2d',Retype='2d',figwidth=24,Ny=100)
        aa.plot_XZ(re_in[-1],'U',0.0,intype='umean_2d',Retype='2d',figwidth=10)
        aa.plot_YZ_xavg(re_in[-1],'V',intype='umean_2d',Retype='2d',figwidth=24,Ny=100)
        aa.plot_YZ_xavg(re_in[-1],'W',intype='umean_2d',Retype='2d',figwidth=24,Ny=100)
# #        aa.plot_YZ_xavg(re_in[-1],'U',intype='SpecReForces',Retype='2d',figwidth=24,Ny=100)
# #        aa.plot_YZ_xavg(re_in[-1],'V',intype='SpecReForces',Retype='2d',figwidth=24,Ny=100)
# #        aa.plot_YZ_xavg(re_in[-1],'W',intype='SpecReForces',Retype='2d',figwidth=24,Ny=100)
        time,Re,Lx,Lz,specF = aa.load_spec(re_in[-1],intype = 'SpecReForces',Retype='2d')
        F,F_ic = aa.div_comp(Lx,Lz,specF)       
        figwidth=24
        aspect=1
        for field in ['U','V','W']:  
            time,Re,Lx,Lz,Y,Z,r = aa.GridXavg(time,Re,Lx,Lz,F_ic,field,Ny=100)
          
            rbar = np.max([-np.min(r),np.max(r)])
          
            plt.figure(figsize=(figwidth,figwidth/3))
            plt.title("F_ic, Re = %.2e, Lx = %.2f, Lz = %.2f, \n Field: %s, time= %.4f" % (Re,Lx,Lz,field,time))
            plt.pcolor(Z,Y,r,vmin=-rbar,vmax=rbar,cmap='bwr')
            ax=plt.gca()
            ax.set_aspect(aspect)
            plt.colorbar()
            plt.ylabel(r'$y$')
            plt.xlabel(r'$z$')
            plt.show()
#
# # figwidths=[10,10,10]
# for ii,run in enumerate(runs):
#     print(run)
#     turb_outs = sorted(glob.glob(idir+run+'/turb*'))
#     for turb in turb_outs[-1:]:
#         aa.plot_turb(turb,figwidth=10)
        
###############
### Spectra ###
###############

#for ii,run in enumerate(runs):
##     print(run)
#    if '_1024' in run:
#        Nx = 1024
#        Nz = 512
#    else:
#        Nx = 512
#        Nz = 256
#    vel_specs = sorted(glob.glob(idir+run+'/vel_sp*'))
#    ks,spec=np.genfromtxt(vel_specs[-1]).T
#    kxs = ks[:int(Nx/2)]
#    kzs = ks[int(Nx/2):]
#    KEx = spec[:int(Nx/2)]
#    KEz = spec[int(Nx/2):]
#    plt.loglog(kxs,KEx,label='KEx')
#    plt.loglog(kzs,KEz,label='KEz')
#    plt.title(run)
#    plt.legend()
#    plt.ylim(np.max([1e-13,np.min([np.min(KEx),np.min(KEz)])]),5*np.max([np.max(KEx),np.max(KEz)]))
#    plt.show()
