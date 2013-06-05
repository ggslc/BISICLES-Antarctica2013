HDF5_HOME=/data/ggslc/opt/BISICLES/hdf5/serial/
NETCDF_HOME=/data/ggslc/opt/BISICLES/netcdf/serial/

HDF5_INC_DIR=$(HDF5_HOME)/include
HDF5_LIB_DIR=$(HDF5_HOME)/lib
HDF5_LIBS = -lhdf5 -lz

NETCDF_INCDIR=$(NETCDF_HOME)/include
NETCDF_LIB_DIR=$(NETCDF_HOME)/lib
NETCDF_LIBS = -L$(NETCDF_LIB_DIR) -L$(HDF5_LIB_DIR) -lnetcdff -lnetcdf -lhdf5_la $(HDF5LIBS)


src:
	cd src && $(MAKE)


