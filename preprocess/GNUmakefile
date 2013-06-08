include Make.defs

VPATH = src:downloadeddata:basedata


FCFLAGS +=  -I. -I $(NETCDF_INCDIR)
LIBFLAGS = $(NETCDF_LIB_FLAGS)

COMOBJS = ncio.o mgrelax.o
MAINS=subsetbm2.f90 prepcontrol.f90 rignot.f90 albtobed.f90
MAINOBJS = $(patsubst %.f90,%.o,$(MAINS))

EXES = $(patsubst %.f90,%.ex,$(MAINS))


all-exes : $(COMOBJS) $(MAINOBJS) $(EXES)

%.ex :  %.o $(OBJS)
	$(FC) $(FCFLAGS)   $< -o $@ $(COMOBJS) $(LIBFLAGS)

%.o : %.f90 
	$(FC) $(FCFLAGS) -c  $< -o $@

all-data : ase-prep antarctica-prep fris-prep

basedata/Antarctica-velocity-1km.nc: rignot.ex src/rawdata.nml downloadeddata/Antarctica_ice_velocity.nc
	./rignot.ex src/rawdata.nml

basedata/Antarctica-bedmap2-1km.nc: subsetbm2.ex src/rawdata.nml downloadeddata/bedmap2_bin.zip 
	./subsetbm2.ex src/rawdata.nml

basedata/Antarctica-temperature-4km.nc basedata/Antarctica-sectors-1km.nc : albtobed.ex src/rawdata.nml downloadeddata/ALBMAP_i2s_4BISICLES.nc
	./albtobed.ex src/rawdata.nml

BASEDATA=basedata/Antarctica-velocity-1km.nc basedata/Antarctica-bedmap2-1km.nc basedata/Antarctica-temperature-4km.nc basedata/Antarctica-sectors-1km.nc


%-geometry-1km.nc %-ctrldata-1km.nc %-temperature-4km.nc %-sectors-1km.nc : src/%.nml $(EXES) $(BASEDATA)
	./prepcontrol.ex $<

%-geometry-1km.2d.hdf5: %-geometry-1km.nc
	$(NCTOAMR) $< $@ thk topg

%-temperature-4km.2d.hdf5: %-temperature-4km.nc
	$(NCTOAMR) $< $@ temp000000 temp000001 temp000002 temp000003 temp000004 temp000005 temp000006 temp000007 temp000008 temp000009

%-ctrldata-1km.2d.hdf5: %-ctrldata-1km.nc
	$(NCTOAMR) $< $@ umod umodc btrc

%-sectors-1km.2d.hdf5: %-sectors-1km.nc
	$(NCTOAMR) $< $@ smask

%-prep: %-geometry-1km.2d.hdf5 %-ctrldata-1km.2d.hdf5 %-sectors-1km.2d.hdf5 %-temperature-4km.2d.hdf5 
	ls -l $? > $@

.PHONY: clean 

clean :
	$(RM) $(COMOBJS) $(MAINOBJS) $(EXES) *.mod


#having gone to the trouble of forging the ring, let us not destroy it
.PRECIOUS: %-geometry-1km.nc %-ctrldata-1km.nc %-temperature-4km.nc %-sectors-1km.nc %-geometry-1km.2d.hdf5 %-ctrldata-1km.2d.hdf5 %-temperature-4km.2d.hdf5 %-sectors-1km.2d.hdf5


