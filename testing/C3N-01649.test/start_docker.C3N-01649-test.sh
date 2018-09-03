# Basic script to start bash in a SomaticWrapper container in standard Docker 

# SomaticWrapper.workflow starts SomaticWrapper docker container with more features
# and works with MGI's non-standard Docker environment

# This starts mwyczalkowski/somatic-wrapper:latest and maps directories:
# Container: /data  
# Host: /Users/mwyczalk/src/SomaticWrapper/data

# Note that DATAD is mapped to the container
# SomaticWrapper work directory is $DATAD/data
#   This allows directories not executed by SomaticWrapper (e.g., A_Reference) to exist on the data partition too
#DATAD="/Users/mwyczalk/Projects/Rabix/TinDaisy/StrelkaDemo.dat"

# Testing using SomaticWrapper results
DATAD="/diskmnt/Projects/Users/hsun/beta_tinDaisy/compare/mgi_sw_C3N-01649"

# results of past TinDaisy run
DATAD="/diskmnt/Projects/Users/hsun/beta_tinDaisy/tmp/tin-daisy.0814/results/TinDaisy.workflow-2018-08-14-210442.362"
IMAGED="/home/mwyczalk_test/data/docker/data"  # a second volume to mount with per-image dagta
GDCDAT="/diskmnt/Projects/cptac_downloads/data/GDC_import" # where GDC BAMs live.  Mounted to /GDC_import

#IMAGE="mwyczalkowski/somatic-wrapper:latest"
IMAGE="cgc-images.sbgenomics.com/m_wyczalkowski/somatic-wrapper:cwl"

docker run -v $DATAD:/data -v $IMAGED:/image -v $GDCDAT:/GDC_import -it $IMAGE

# To start another terminal in running container, first get name of running container with `docker ps`,
# then start bash in it with,
# `docker exec -it <container_name> bash`

