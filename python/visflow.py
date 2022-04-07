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


##########################
### Flow visualization ###
##########################

# figwidths=[10,10,10]
for ii,run in enumerate(runs):
    print(run)
    state_outs = sorted(glob.glob(idir+run+'/state*'))
    print(state_outs[:])
    state = input("Which state number to plot?")
    field = input("Which field? U,V,W ")
    yloc = input("Which ylocation? ")
    yloc = float(yloc)
    state = str(state).zfill(4)
    aa.plot_state(idir+run+'/state'+state+'.cdf.dat',field,yloc,figwidth=10)
        
