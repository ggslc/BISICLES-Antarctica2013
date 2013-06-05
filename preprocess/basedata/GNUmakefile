

all: Antarctica-velocity-1km.nc Antarctica-temperature-4km.nc Antarctica-bedmap2-1km.nc

Antarctica-velocity-1km.nc: ../src/rignot.ex ../downloadeddata/Antarctica_ice_velocity.nc
	../src/rignot.ex ../src/rawdata.nml

Antarctica-temperature-4km.nc: ../src/albtobed.ex ../downloadeddata/ALBMAP_i2s_4BISICLES.nc
	../src/albtobed.ex ../src/rawdata.nml

Antarctica-bedmap2-1km.nc: ../src/subsetbm2.ex ../downloadeddata/bedmap2_bin.zip
	../src/subsetbm2.ex ../src/rawdata.nml


