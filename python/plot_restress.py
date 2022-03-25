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
    plt.semilogy(t,ke,label=run)
    
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
#     print(state_outs[:])
    for state_out in np.array(state_outs)[:1]:
        print(state_out)
        aa.plot_state(state_out,'U',0.0,figwidth=10)

# figwidths=[10,10,10]
for ii,run in enumerate(runs):
    print(run)
    state_outs = sorted(glob.glob(idir+run+'/state*'))
#     print(state_outs[:])
    for state_out in np.array(state_outs)[-1:]:
        print(state_out)
        aa.plot_state(state_out,'U',0.0,figwidth=10)

# Even modes
field = input("Which (even) field to look at? (options: umean, wmean, uumean, uwmean, wwmean, vvmean)")
for ii,run in enumerate(runs):
    print(run)
    state_outs = sorted(glob.glob(idir+run+'/restress_even*'))
#     print(state_outs[:])
    for state_out in np.array(state_outs)[-1:]:
#         print(state_out)
        aa.plot_restress(state_out,field,Ny=50)

# Odd modes
field = input("Which (odd) field to look at? (options: vmean, uvmean, wvmean)")
for ii,run in enumerate(runs):
    print(run)
    state_outs = sorted(glob.glob(idir+run+'/restress_odd*'))
#     print(state_outs[:])
    for state_out in np.array(state_outs)[-1:]:
#         print(state_out)
        aa.plot_restress(state_out,field,Ny=50)

# Forces modes
direction = input("Plotting forces now. Which force direction? x,y,z? ")
for ii,run in enumerate(runs):
    print(run)
    state_outs = sorted(glob.glob(idir+run+'/restress_odd*'))
    print(state_outs[:])
    outnum = input("Which output? ")
    aa.plot_reforces(idir+run,outnum,direction,Ny=50)
