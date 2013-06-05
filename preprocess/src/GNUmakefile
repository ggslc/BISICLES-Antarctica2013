FC = gfortran
CC=gcc
UNAMEN := $(shell uname -n)

#netcdf location
NETCDF_HOME = /data/ggslc/opt/BISICLES/netcdf/serial
HDF5_HOME =  /data/ggslc/opt/BISICLES/hdf5/serial


FFLAGS = -g  -I. -I$(NETCDF_HOME)/include -I$(HDF5_HOME)/include 
LIBFLAGS = -L$(NETCDF_HOME)/lib -lnetcdff -lnetcdf -L$(HDF5_HOME)/lib -lhdf5_hl -lhdf5 -lz


OBJS = ncio.o mgrelax.o prepcontrol.o

prepcontrol: $(OBJS)
	$(FC) $(FFLAGS) -o prepcontrol $(OBJS)  $(LIBFLAGS)

%.o : %.f90 
	$(FC) $(FFLAGS) -c  $< -o $@

clean: 
	$(RM) $(OBJS) prepcontrol

realclean:
	$(RM) $(OBJS) prepcontrol *.mod
