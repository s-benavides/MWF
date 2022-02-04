#!/bin/bash

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_mit_raffaele    (32)    2-10-1  12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00

# Job
#SBATCH --partition=sched_mit_hill
##SBATCH --partition=sched_any_quicktest
##SBATCH --partition=newnodes
#SBATCH --nodes=2                 ## 256^2 = 8 cores, 512^2 = 16 cores
#SBATCH --ntasks-per-node=12
##SBATCH --mem-per-cpu=3000
##SBATCH --time=0:15:00
#SBATCH --time=12:00:00
#SBATCH -J Lx180Lz80  # sensible name for the job

## load up the correct modules, if required
. /etc/profile.d/modules.sh
module load engaging/openmpi/2.0.3

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Run scripts
mpiexec -n 24 -mca btl_tcp_if_include ib0 ./randIC.out
mv state0000.cdf.dat state.cdf.in
mpiexec -n 24 -mca btl_tcp_if_include ib0 ./main.out



