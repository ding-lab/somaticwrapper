# Index the strelka demo reference
# Note that the file /data/A_Reference/demo20.fa should have already been created
# and modified to use standard reference names (somaticwrapper/docker/2_setup_data.sh does this)


OUTD="/data/A_Reference"
F="demo20.fa"

if [ ! -f $OUTD/$F ]; then
    echo "Error: Reference $OUTD/$F not found.  Please create it using somaticwrapper/docker/2_setup_data.sh"
    exit
fi

echo Preparing reference $F
bash ./prepare_reference.sh $OUTD/$F

