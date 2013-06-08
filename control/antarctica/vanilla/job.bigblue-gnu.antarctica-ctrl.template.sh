#!/bin/bash
#
# job script template for bluecrystal phase 2
#
#PBS -l walltime=@HOURS:00:00,nodes=@NODES:ppn=8
#request  nodes ,ppn processor per node
#
#PBS -j oe
# 
#get the job id no
export JOBNO="`echo $PBS_JOBID | sed s/.master.ic.cluster//`"

#create job director
export BASEDIR=/gpfs/cluster/geog/ggslc/Antarctica2013/Antarctica/control-vanilla/
export DIR=$BASEDIR/@MAXLEVlev
mkdir -p $DIR
cd $DIR

export BASECFG=inputs.@NAME
INFILE=$BASECFG.$JOBNO
echo "INFILE = $INFILE"
cp $BASEDIR/scripts/$BASECFG $INFILE

#ensure that the job exits on checkpoint
#echo "amr.check_exit = true" >> $INFILE


#work out what the latest checkpoint file is (if it exists)
if test -n "$(find . -maxdepth 1 -name 'chk*' -print -quit)"
    then
    LCHK=`ls -th chk* | head -n 1`
    echo "" >> $INFILE #ensure a line break
    echo "amr.restart_file=$LCHK" >> $INFILE
fi
export MYEXPATH="/gpfs/cluster/geog/ggslc/BISICLES/BISICLES/code/controlproblem"
export MYEXEX="control2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex"
export MYEX=$MYEXPATH/$MYEXEX
export CMD="$MYEX $INFILE"
echo "CMD = $CMD"
#
echo "PBS_NODEFILE = $PBS_NODEFILE"
nodelist=`cat $PBS_NODEFILE`



#
echo "JOBNO = $JOBNO" 
echo "nodelist = $nodelist"
#
#   name the configuration file
export CONFILE="$DIR/ib.$JOBNO.conf"
echo "CONFILE =  $CONFILE"
#
#   put nodes in configuration file
#
#number of cores on each node
export cpn=8
#
for iii in `cat $PBS_NODEFILE | uniq` 
do 
  j=0
  while [ $j -lt $cpn ]
  do
        echo $iii >> $CONFILE
	let "j=j+1"
  done
  ssh $iii mkdir -p /local/ggslc
  ssh $iii killall -s KILL $MYEXEX
done



#
#   get the number of processors
NUMPROC=`cat $CONFILE | wc -l`
export NUMPROC
echo $NUMPROC processors
#
#    run job
#
if mpirun -np $NUMPROC -machinefile $CONFILE -x CH_TIMER=1  $CMD;
then
    for iii in `cat $PBS_NODEFILE | uniq`;do 
	ssh $iii "zip $MYDIR/$JOBNO.zip /local/ggslc/pout.* && rm /local/ggslc/pout.*" 
    done
    zip $JOBNO.zip time.table.*
else
    exit 1
fi
