#!/bin/bash

# run a given SomaticWrapper step (or set of steps) via LSF scheduler on MGI
# Usage: submit_SomaticWrapper_MGI.sh [options] STEP CONFIG 
#   STEP is the step number in SomaticWrapper.pl.  Names (e.g., "parse_pindel") can also be used (untested)
#   CONFIG is the configuration file as described in SomaticWrapper.pl, with a path that is visible
#       from within the running container (e.g., /data/foo)

# Optional arguments:
# -m MEMGb - integer indicating number of gigabytes to allocate.  Default value (set by MGI) is possibly 8
# -D DATAD - path to container's /data is mounted on host [default /gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data]
# -h DOCKERHOST - define a host to execute the image
# -s SCRIPTD_BASE - Script run base directory, where bsub output will be written.  $SAMPLE_NAME is appended
#    [default: /gscuser/mwyczalk/projects/SomaticWrapper/runtime_bsub/]
# -S SAMPLE_NAME - sample name, here used only to modify SCRIPTD path.  Required
# -p SW - Path to SomaticWrapper base dir  [default: /gscuser/mwyczalk/projects/SomaticWrapper/somaticwrapper]
# -d - dry run, print out run command but don't execute

if [ "$#" -lt 2 ]
then
    echo "Error - invalid number of arguments"
    echo "Usage: $0 [options] STEP CONFIG "
    exit 1
fi

# Set defaults here
# Where container's /data is mounted on host
DATAD="/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data"
BSUB="bsub"
SW="/gscuser/mwyczalk/projects/SomaticWrapper/somaticwrapper"
SCRIPTD_BASE="/gscuser/mwyczalk/projects/SomaticWrapper/runtime_bsub"

# http://wiki.bash-hackers.org/howto/getopts_tutorial
while getopts ":m:D:h:s:S:p:d" opt; do
  case $opt in
    m)
      MEMGB=$OPTARG  
      echo "Setting memory $MEMGB Gb" >&2
      # MEM will be undefined (or blank) if not set here
      MEM="-R \"rusage[mem=${MEMGB}000]\" -M ${MEMGB}000000"
      ;;
    D)
      DATAD=$OPTARG
      ;;
    h)
    # User may define a specific host to avoid downloading image each run
    # This option will be undefined if not set here
      DOCKERHOST="-m $OPTARG"
      echo "Setting host $OPTARG" >&2
      ;;
    s)
      SCRIPTD_BASE=$OPTARG
      echo "Setting script base directory $SCRIPTD_BASE" >&2
      ;;
    S)      # this is required
      SAMPLE_NAME=$OPTARG
      ;;
    P)
      SW=$OPTARG
      ;;
    d)
      echo "Dry run" >&2
      BSUB="echo bsub"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z $SAMPLE_NAME ]; then
echo Sample Name not defined \(-S SAMPLE_NAME\)
exit 1
fi

# /gscuser/mwyczalk/projects/SomaticWrapper/SW_testing/BRCA77/1_start_runs.sh

shift $((OPTIND-1))
STEP=$1; shift
echo Running step $STEP >&2

# Config file will be read from within container, so the path must be visible from it
# One complication is that we cannot parse $CONFIG from host (/data not mapped); 
#CONFIG=$2
CONFIG=$1; shift
echo Using configuration file $CONFIG >&2

echo "/data mounts to $DATAD " >&2
export LSF_DOCKER_VOLUMES="$DATAD:/data"

# Here we're putting together a script which will be run in new container to launch a job after sourcing environment variables
SCRIPTD="$SCRIPTD_BASE/$SAMPLE_NAME"
mkdir -p $SCRIPTD
echo bsub run output directory is $SCRIPTD >&2
SCRIPT="$SCRIPTD/bsub_run-step_$STEP.sh"
echo Generating run script $SCRIPT >&2

cat << EOF > $SCRIPT
#!/bin/bash

# This is an automatically generated script for launching bsub jobs

source /home/bps/mgi-bps.bashrc
cd $SW
perl $SW/SomaticWrapper.pl /data/data $STEP $CONFIG
EOF

# logs will be written to $SCRIPTD/bsub_run-step_$STEP.err, .out
LOGS="-e $SCRIPTD/bsub_run-step_$STEP.err -o $SCRIPTD/bsub_run-step_$STEP.out"
rm -f $SCRIPTD/bsub_run-step_$STEP.err $SCRIPTD/bsub_run-step_$STEP.out


$BSUB -q research-hpc $DOCKERHOST $LOGS $MEM -a 'docker (mwyczalkowski/somatic-wrapper:mgi)'  "/bin/bash $SCRIPT"
