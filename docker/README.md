Workflow in this directory is created for development and testing of the docker image.  All
scripts run on host.  Testing is done with the Strelka test dataset

Also note that the Dockerfile has passwords associated with it, so it is not for distribution until 
somatic_wrapper is public

TODO: Split the development stuff (Docker, etc) from things which need to be run to initialize the
image on the host.
