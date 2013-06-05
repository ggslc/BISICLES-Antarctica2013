include ../Make.defs

FCFLAGS +=  -I. -I $(NETCDF_INCDIR)
LIBFLAGS = $(NETCDF_LIB_FLAGS)

COMOBJS = ncio.o mgrelax.o
MAINS=subsetbm2.f90 prepcontrol.f90 rignot.f90 albtobed.f90
MAINOBJS = $(patsubst %.f90,%.o,$(MAINS))

EXES = $(patsubst %.f90,%.ex,$(MAINS))


all : $(COMOBJS) $(MAINOBJS) $(EXES)

%.ex :  %.o $(OBJS)
	$(FC) $(FCFLAGS)   $< -o $@ $(COMOBJS) $(LIBFLAGS)

%.o : %.f90 
	$(FC) $(FCFLAGS) -c  $< -o $@

clean :
	$(RM) $(COMOBJS) $(MAINOBJS) $(EXES) *.mod

