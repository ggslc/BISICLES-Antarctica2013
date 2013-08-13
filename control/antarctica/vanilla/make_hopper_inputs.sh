

mkin()
{
    NAME=antarctica-ctrl-"$lev"lev    
    SCATCHDIR="\/scratch\/scratchdirs\/cornford\/bigger_pigthwaites"
    outdir="$pthwdir\/$subdir\/$NAME"
    SUBS="-e s/@SCRATCHDIR/$SCRATCHDIR/ -e s/@NAME/$NAME/ -e s/@QUEUE/$QUEUE/ -e s/@HOURS/$HOURS/ -e s/@MPPWIDTH/$MPPWIDTH/ -e s/@LEV/$LEV/"
    sed $SUBS inputs.template > scripts/inputs.$NAME
    sed $SUBS job.hopper.template.sh > scripts/job.hopper.$NAME.sh
   
}


LEV=5
HOURS=1
QUEUE=regular
MPPWIDTH=144
mkin