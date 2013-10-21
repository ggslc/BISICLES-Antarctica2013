FLATTEN=/data/ggslc/opt/BISICLES/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.OPT.ex
EXTRACT=/data/ggslc/opt/BISICLES/BISICLES/code/filetools/extract2d.Linux.64.g++.gfortran.DEBUG.OPT.ex

exflat()
{
   echo $CONTROL
   $EXTRACT $DIR/$FILE tmp.2d.hdf5 $EXARGS
   $FLATTEN tmp.2d.hdf5 ase-data-$CTRL-1km.2d.hdf5 2 
   rm tmp.2d.hdf5
}

DIR=../control/ase/mass-balance/test5/
FILE=ase-test5-ctrl-2lev-outer.000128.2d.hdf5
CTRL=mb
EXARGS="thck topg Cwshelf C muCoef divuho" 

exflat


