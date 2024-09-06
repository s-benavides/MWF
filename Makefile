INSTDIR		= ./install/
PROGDIR		= ./program/
UTILDIR		= ./utils/
UTIL		= randIC

TRANSFORM	= fftw5
MODES		= M4
MODSOBJ		= mpi.o parameters.o \
		modes.o variables.o transform.o velocity.o \
		turb.o io.o 

NETCDFPATH = /media/apps/avx512-2021/software/netCDF-Fortran/4.6.1-gompi-2023a
FFTWPATH = /media/apps/avx512-2021/software/FFTW.MPI/3.3.10-gompi-2023a
MPIPATH = /media/apps/avx512-2021/software/OpenMPI/4.1.5-GCC-12.3.0

COMPILER        = /media/apps/avx512-2021/software/GCCcore/13.2.0/bin/gfortran
COMPFLAGS       = -fallow-argument-mismatch -cpp -traditional-cpp -ffree-line-length-none -x f95-cpp-input -fPIC -O3 -march=native -ffast-math -funroll-loops -c \
		-I$(MPIPATH)/include/ -I$(MPIPATH)/lib/ \
		-I$(NETCDFPATH)/include -I/media/apps/avx512-2021/software/bzip2/1.0.8-GCCcore-12.3.0/include -I/media/apps/avx512-2021/software/netCDF/4.9.2-gompi-2023a/include -DgFortran \
		-I$(FFTWPATH)/include/

LIBS            = $(FFTWPATH)/lib/libfftw3.a -L$(FFTWPATH)/lib/ \
		$(NETCDFPATH)/lib/libnetcdff.a -L$(NETCDFPATH)/lib/ -L$(NETCDFPATH)/lib64/ -lnetcdf \
		-L$(MPIPATH)/lib -L/media/apps/avx512-2021/software/hwloc/2.9.1-GCCcore-12.3.0/lib -L/media/apps/avx512-2021/software/libevent/2.1.12-GCCcore-12.3.0/lib -Wl,-rpath -Wl,$(MPIPATH)/lib -Wl,-rpath -Wl,/media/apps/avx512-2021/software/hwloc/2.9.1-GCCcore-12.3.0/lib -Wl,-rpath -Wl,/media/apps/avx512-2021/software/libevent/2.1.12-GCCcore-12.3.0/lib -Wl,--enable-new-dtags -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi  ## 'mpif90 -show' output

#------------------------------------------------------------------------

all : 	$(MODSOBJ) $(PROGDIR)main.f90
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)main.f90
	$(COMPILER) -o ./main.out main.o $(MODSOBJ) $(LIBS)

install : main.out
	if test ! -d $(INSTDIR); then mkdir -p $(INSTDIR); fi
	mv ./main.out $(INSTDIR)
	date > $(INSTDIR)/main.info
	echo $(HOSTNAME) >> $(INSTDIR)/main.info
	pwd >> $(INSTDIR)/main.info
	echo $(COMPILER) $(COMPFLAGS) >> $(INSTDIR)/main.info
	echo $(MODES) >> $(INSTDIR)/main.info
	grep define parallel.h | grep _Np >> $(INSTDIR)/main.info
	cut -d! -f1 $(PROGDIR)parameters.f90 | grep = | \
	   cut -d: -f3  >> $(INSTDIR)/main.info

util : 	$(MODSOBJ) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) $(COMPFLAGS) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) -o ./$(UTIL).out $(UTIL).o $(MODSOBJ) $(LIBS)

#------------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.il core *.out

#------------------------------------------------------------------------
netcdf_int.o : $(PROGDIR)netcdf_int.f90 
	$(COMPILER) $(COMPFLAGS) -o netcdf_int.o $(PROGDIR)netcdf_int.f90

io.o : $(PROGDIR)io.f90 
	$(COMPILER) $(COMPFLAGS) -o io.o $(PROGDIR)io.f90

turb.o : $(PROGDIR)turb.f90 velocity.o mpi.o 
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)turb.f90

mpi.o : $(PROGDIR)mpi.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)mpi.f90

parameters.o :	$(PROGDIR)parameters.f90 parallel.h
		$(COMPILER) $(COMPFLAGS) $(PROGDIR)parameters.f90

timestep.o : $(PROGDIR)timestep.f90 variables.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)timestep.f90

transform.o : $(PROGDIR)transform.$(TRANSFORM).f90 variables.o
	$(COMPILER) $(COMPFLAGS) -o transform.o \
	$(PROGDIR)transform.$(TRANSFORM).f90

variables.o : $(PROGDIR)variables.f90 parameters.o mpi.o
	$(COMPILER) $(COMPFLAGS) -o variables.o $(PROGDIR)variables.f90

velocity.o : $(PROGDIR)velocity.f90 transform.o modes.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)velocity.f90

modes.o : $(PROGDIR)modes.$(MODES).f90 parameters.o
	$(COMPILER) $(COMPFLAGS) -o modes.o \
	$(PROGDIR)modes.$(MODES).f90
