
DOMAIN=ase
BASEDIR=`pwd`/../../ 
PREPDIR=$BASEDIR/preprocess
SHAREDIR=$BASEDIR/share
RUNPARENT=`pwd`

mkin()
{
    tc=$(( lev - 1 ))
    NAME=ase-spinup.$lev"lev"
    OUTDIR="$RUNPARENT/$NAME"
    echo $OUTDIR
    mkdir -p $OUTDIR
    SUBS="-e s:@PREPDIR:$PREPDIR: -e s:@SHAREDIR:$SHAREDIR: -e s:@OUTDIR:$OUTDIR: -e s:@NAME:$NAME: -e s:@TAGCAP:$tc: -e s:@LEV:$lev:" 
    sed $SUBS inputs.template > scripts/inputs.$NAME 
}

mkdir -p scripts

for lev in 0 1 2 3 4; do
    mkin
done



