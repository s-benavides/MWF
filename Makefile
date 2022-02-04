INSTDIR		= ./install/
PROGDIR		= ./program/
UTILDIR		= ./utils/
UTIL		= randIC

TRANSFORM	= fftw5
MODES		= M4
MODSOBJ		= mpi.o parameters.o \
		modes.o variables.o transform.o velocity.o \
		turb.o io.o 

NETCDFPATH = /home/santiago_b/.local
##NETCDFPATH = /home/software/netcdf/4.6.3
FFTWPATH = /home/santiago_b/FFTW3_GHOST
MPIPATH = /cm/shared/engaging/openmpi/2.0.3

COMPILER        = mpif90
##COMPFLAGS       = -ffree-line-length-none -x f95-cpp-input -c -O3 \

COMPFLAGS       = -cpp -traditional-cpp -O3 -funroll-loops -ffree-line-length-none -x f95-cpp-input -c \
		-I$(MPIPATH)/include/ -I$(MPIPATH)/lib/\
		-I$(NETCDFPATH)/include/ \
		-I$(FFTWPATH)/include/
LIBS            = -lm \
		-L$(NETCDFPATH)/lib -lnetcdff -L$(NETCDFPATH)/lib -lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl \
		-L$(FFTWPATH)/lib -lfftw3 \
		-L$(MPIPATH)/lib/ \
		-pthread -Wl,-rpath -Wl,$(MPIPATH)/lib -Wl,--enable-new-dtags -lmpi_usempi -lmpi_mpifh -lmpi

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
