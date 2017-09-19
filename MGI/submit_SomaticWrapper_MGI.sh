#!/bin/bash

if [ "$#" -ne 2 ]
then
    echo "Error - invalid number of arguments"
    echo "Usage: $0 STEP CONFIG"
    exit 1
fi

STEP=$1
echo Running step $STEP

# sw.config is the Strelka Demo
#CONFIG="$DATAD/data/SWtest/sw.config"
# Config file will be read from within container, so the path must be visible from it
CONFIG=$2


# Where container's /data is mounted on host
DATAD="/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data"

echo Mapping /data to $DATAD

export LSF_DOCKER_VOLUMES="$DATAD:/data"

# Root of somaticwrapper project on host
SW="/gscuser/mwyczalk/projects/SomaticWrapper/somaticwrapper"

# Stick to a host where we've already downloaded the image for performance reasons
HOST="-m blade18-1-1.gsc.wustl.edu"

# Here we're putting together a script which will be run in new container to launch a job after sourcing environment variables
SCRIPTD="/gscuser/mwyczalk/projects/SomaticWrapper/runtime_bsub"
mkdir -p $SCRIPTD
SCRIPT="$SCRIPTD/bsub_run-step_$STEP.sh"

cat << EOF > $SCRIPT
#!/bin/bash

# This is an automatically generated script for launching bsub jobs

source /home/bps/mgi-bps.bashrc
cd $SW
perl $SW/SomaticWrapper.pl /data/data $STEP $CONFIG
EOF

echo Written configuration to $SCRIPT

# logs will be written to $SCRIPTD/bsub_run-step_$STEP.err, .out
LOGS="-e $SCRIPTD/bsub_run-step_$STEP.err -o $SCRIPTD/bsub_run-step_$STEP.out"


# To run docker container in background, do something like,
bsub -q research-hpc $HOST $LOGS -a 'docker (mwyczalkowski/somatic-wrapper:mgi)'  "/bin/bash $SCRIPT"
