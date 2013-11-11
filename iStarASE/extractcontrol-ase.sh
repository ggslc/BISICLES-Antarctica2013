FLATTEN=/data/ggslc/opt/BISICLES/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.OPT.ex
EXTRACT=/data/ggslc/opt/BISICLES/BISICLES/code/filetools/extract2d.Linux.64.g++.gfortran.DEBUG.OPT.ex
PYTHONF=/data/ggslc/opt/BISICLES/BISICLES/code/filetools/pythonf2d.Linux.64.g++.gfortran.DEBUG.OPT.ex

exflat()
{
   echo $CONTROL
   $EXTRACT $DIR/$FILE tmp.2d.hdf5 $EXARGS
   $FLATTEN tmp.2d.hdf5 tmp2.hdf5 2
   if [ $CTRL == mb ] 
   then
       
       TUPLE="thck,topg,Cwshelf,C,muCoef,divuho"
       $PYTHONF tmp2.hdf5 ase-data-$CTRL-1km.2d.hdf5 asepostcontrol asepostcontrol  $TUPLE $TUPLE
   else
       cp tmp2.hdf5 ase-data-$CTRL-1km.2d.hdf5
   fi
   rm tmp.2d.hdf5 tmp2.hdf5
}

export PYTHONPATH=`pwd`

DIR=../control/ase/mass-balance/test7/
FILE=ase-test7-ctrl-2lev-outer.000512.2d.hdf5
CTRL=mb
EXARGS="thck topg Cwshelf C muCoef divuho" 
exflat

DIR=../control/ase/mass-balance/test7/
FILE=ase-test7-ctrl-2lev-outer.000000.2d.hdf5
CTRL=zero
EXARGS="thck topg Cwshelf C muCoef divuho" 


exflat


