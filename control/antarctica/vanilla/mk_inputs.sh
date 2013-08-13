

getnodes()
{
    case $MAXLEV in
	0)
	    NODES=1; HOURS=8;;
	1)
	    NODES=1; HOURS=8;;
	2)
	    NODES=1; HOURS=8;;
	3)
	    NODES=2; HOURS=8;;
	4)
	    NODES=4; HOURS=8;;
    esac


}

EXPT=antarctica-ctrl

mkdir -p scripts

for MAXLEV in 2 3 4; do

    getnodes;

    for CONTROL in vanilla ; do
	NAME=$EXPT"-"$MAXLEV"lev"
	#PWD=`pwd`
	SUBS="-e s/@EXPT/$EXPT/ -e s/@LEV/$MAXLEV/ -e s/@NAME/$NAME/ -e s:@PWD:$PWD: "
	sed $SUBS inputs.$EXPT.template > scripts/inputs.$NAME	
	sed $SUBS -e s/@NODES/$NODES/ -e s/@HOURS/$HOURS/ job.bigblue-gnu.$EXPT.template.sh > scripts/job.bigblue-gnu.$NAME.sh
	sed $SUBS -e s/@NODES/$NODES/ -e s/@HOURS/$HOURS/ job.bigblue-intel.$EXPT.template.sh > scripts/job.bigblue-intel.$NAME.sh
    done
done

EXPT=antarctica-petsc-ctrl

for MAXLEV in 4; do

    getnodes;

    for CONTROL in vanilla ; do
	NAME=$EXPT"-"$MAXLEV"lev"
	#PWD=`pwd`
	SUBS="-e s/@EXPT/$EXPT/ -e s/@MAXLEV/$MAXLEV/ -e s/@NAME/$NAME/ -e s:@PWD:$PWD: "
	sed $SUBS inputs.$EXPT.template > scripts/inputs.$NAME
	sed $SUBS -e s/@NODES/$NODES/ -e s/@HOURS/$HOURS/ job.bigblue-gnu.$EXPT.template.sh > scripts/job.bigblue-gnu.$NAME.sh
	
    done
done