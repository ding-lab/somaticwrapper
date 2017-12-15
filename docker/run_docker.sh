# Basic script to start bash in a SomaticWrapper container in standard Docker 

# SomaticWrapper.workflow starts SomaticWrapper docker container with more features
# and works with MGI's non-standard Docker environment

# This starts mwyczalkowski/somatic-wrapper:latest and maps directories:
# Container: /data  
# Host: /Users/mwyczalk/src/SomaticWrapper/data

# Note that DATAD is mapped to the container
# SomaticWrapper work directory is $DATAD/data
#   This allows directories not executed by SomaticWrapper (e.g., A_Reference) to exist on the data partition too
DATAD="$HOME/src/SomaticWrapper/data"
IMAGE="mwyczalkowski/somatic-wrapper:latest"

docker run -v $DATAD:/data -it $IMAGE

# To start another terminal in running container, first get name of running container with `docker ps`,
# then start bash in it with,
# `docker exec -it <container_name> bash`

