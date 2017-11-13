#!/bin/bash

if [ "$#" -ne 1 ]
then
    echo "Error - invalid number of arguments"
    echo "Usage: $0 command "
    exit 1
fi

CMD=$1
echo Running $CMD

# Where container's /data is mounted on host
DATAD="/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data"

echo Mapping /data to $DATAD

export LSF_DOCKER_VOLUMES="$DATAD:/data"

# Stick to a host where we've already downloaded the image for performance reasons
#DOCKERHOST="-m blade18-1-1.gsc.wustl.edu"

# logs will be written to ./runtime/bsub_run-DATE.err, .out
NOW=$(date +%Y%m%d_%H%M%S)
OUTD="runtime"
mkdir -p $OUTD
RUND=$(readlink -f $OUTD)  # full path to this directory so visible on container
LOGS="-e $RUND/bsub_run-step_$NOW.err -o $RUND/bsub_run-step_$NOW.out"
echo Writing logs to $LOGS

# Extra memory required for Pindel run
# Ideally, should specify memory requirements in configuration file

#MEMGB=30
#MEM="-R \"rusage[mem=${MEMGB}000]\" -M ${MEMGB}000000"

# Here we're putting together a script which will be run in new container to launch a job after sourcing environment variables
# This script called ./runtime/run_$NOW.sh
SCRIPT="$RUND/run_$NOW.sh"

cat << EOF > $SCRIPT
#!/bin/bash

# This is an automatically generated script for launching bsub jobs

source /home/bps/mgi-bps.bashrc
$CMD 

EOF

echo Written configuration to $SCRIPT

bsub -q research-hpc $DOCKERHOST $LOGS $MEM -a 'docker (mwyczalkowski/somatic-wrapper:mgi)'  "/bin/bash $SCRIPT"
