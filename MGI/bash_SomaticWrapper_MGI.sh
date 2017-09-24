# make sure ./somaticwrapper is up to date

# echo Be sure to do,
# echo source /home/bps/mgi-bps.bashrc

HOST="-m blade18-1-1.gsc.wustl.edu"

# Where container's /data is mounted on host
DATAD="/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data"
echo Mapping /data to $DATAD
export LSF_DOCKER_VOLUMES="$DATAD:/data"


bsub -q research-hpc $HOST -Is -a "docker(mwyczalkowski/somatic-wrapper:mgi)" "/bin/bash /home/bps/mgi-bps_start.sh"
#COPY image-init/mgi-bps_start.sh /home/bps


# CONFIG="/gscuser/mwyczalk/projects/SomaticWrapper/config/01BR001.config"
# To run docker container in background, do something like,
# bsub -q research-hpc -a 'docker (mwyczalkowski/somatic-wrapper:mgi)'  "source /home/bps/mgi-bps.bashrc; perl SomaticWrapper.pl /data/data 2 ../config/01BR001.config "
