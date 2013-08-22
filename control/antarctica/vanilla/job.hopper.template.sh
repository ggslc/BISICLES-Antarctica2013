#!/bin/bash
#
# job script template for hopper
#PBS -q @QUEUE
#PBS -l mppwidth=@MPPWIDTH
#PBS -l walltime=@HOURS:00:00
#PBS -j oe
# 
export SRCDIR=/global/homes/c/cornford/bisicles-antarctica2013/control/antarctica/vanilla/@SUBDIR
export RUNDIR=@SCRATCHDIR/@SUBDIR/@NAME
mkdir -p $RUNDIR
cd $RUNDIR
export BASECFG=inputs.@NAME
INFILE=$BASECFG
cp $SRCDIR/scripts/$BASECFG $INFILE

export EX="/global/homes/c/cornford/BISICLES/BISICLES/code/controlproblem/control2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.ex"
export CMD="$EX $INFILE"
echo "CMD = $CMD"
#run the job
aprun -n @MPPWIDTH $CMD