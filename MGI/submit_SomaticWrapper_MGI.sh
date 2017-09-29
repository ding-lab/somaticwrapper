#!/bin/bash

if [ "$#" -lt 3 ]
then
    echo "Error - invalid number of arguments"
    echo "Usage: $0 STEP CONFIG SAMPLE_NAME [MEMGb]"
    exit 1
fi

# STEP is the step number in SomaticWrapper.pl.  Names (e.g., "parse_pindel") can also be used (untested)
# CONFIG is the configuration file as described in SomaticWrapper.pl, with a path that is visible
#   from within the running container (e.g., /data/foo)
# SAMPLE_NAME is sample_name of this run as described in the CONFIG file.  
# MEMGb is an integer indicating number of gigabytes to allocate.  Default value (set by MGI) is possibly 8

STEP=$1
echo Running step $STEP

# sw.config is the Strelka Demo
#CONFIG="$DATAD/data/SWtest/sw.config"
# Config file will be read from within container, so the path must be visible from it
CONFIG=$2

# One complication is that we cannot parse $CONFIG from host (/data not mapped); SAMPLE_NAME here is used only 
# to define paths of run output files
SAMPLE_NAME=$3

MEMGB=$4
if [ ! -z $MEMGB ]; then
# Extra memory required for Pindel run
MEM="-R \"rusage[mem=${MEMGB}000]\" -M ${MEMGB}000000"
fi

# Where container's /data is mounted on host
DATAD="/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data"

echo Mapping /data to $DATAD

export LSF_DOCKER_VOLUMES="$DATAD:/data"

# Root of somaticwrapper project on host
SW="/gscuser/mwyczalk/projects/SomaticWrapper/somaticwrapper"

## Stick to a host where we've already downloaded the image for performance reasons
#DOCKERHOST="-m blade18-1-1.gsc.wustl.edu"

# Here we're putting together a script which will be run in new container to launch a job after sourcing environment variables
SCRIPTD="/gscuser/mwyczalk/projects/SomaticWrapper/runtime_bsub/$SAMPLE_NAME"
mkdir -p $SCRIPTD
echo bsub run output directory is $SCRIPTD
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
rm -f $SCRIPTD/bsub_run-step_$STEP.err $SCRIPTD/bsub_run-step_$STEP.out


bsub -q research-hpc $DOCKERHOST $LOGS $MEM -a 'docker (mwyczalkowski/somatic-wrapper:mgi)'  "/bin/bash $SCRIPT"
