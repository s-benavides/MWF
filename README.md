This code simulated Model Waleffe flow in a retangular domain.

To build:
1. Edit Makefile to point to netcdf & fftw
2. make;make install

If starting fresh, a random initial condition will be generated. If state relaminarises then increase d_E0 & recompile, or run at higher Re.

To run.
1. Copy main.info, main.out & parameter.inp to folder.
2. Run ./main.out
3. vel_energy.dat keeps a running output of the total energy.
3. Delete RUNNING file to softly kill run.
4. For parallel simulation edit _Np in parallel.h & recompile.

To control.
1. program/parameters.f90 contains information regarding resolution.
2. main.info has a copied of the parameters used at compilation.
3. Recompile.
4. Parameter.inp contains parameters such as Re, Lx, Lz, and more. These don't require you to recompile and can be changed in the local run directory.

To plot output.
1. In matlab [x,z,u]=GridUy('state.cdf.in','U',0.) extracts the U field at y=0.
2. contourf(x,z,u) to visualise.

Questions?
Feel free to contact me.
