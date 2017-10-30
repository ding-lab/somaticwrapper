# make sure ./somaticwrapper is up to date

# echo Be sure to do,
# echo source /home/bps/mgi-bps.bashrc

# Define HOST if want to use the same host each time.  This is useful to keep docker image from being downloaded 
# every time (an individual machine will keep a cached image).  However, intensive processing may overload a
# single machine
#DOCKERHOST="-m blade18-1-1.gsc.wustl.edu"

# Where container's /data is mounted on host
DATAD="/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data"
echo Mapping /data to $DATAD
export LSF_DOCKER_VOLUMES="$DATAD:/data"

MEMGB=24
MEM="-R \"rusage[mem=${MEMGB}000]\" -M ${MEMGB}000000"

# Note memory usage, hoping to keep pindel from dying
bsub -q research-hpc $DOCKERHOST -Is $MEM -a "docker(mwyczalkowski/somatic-wrapper:mgi)" "/bin/bash /home/bps/mgi-bps_start.sh"


