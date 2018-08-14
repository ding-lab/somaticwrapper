Build docker image containing SomaticWrapper (CWL branch) and all
necessary software to run it.

`./docker.testing` includes scripts to start docker image and execute
specific steps from within the docker environment.  This is useful for
debugging and feature development.

`./StrelkaDemo.dat` includes small test dataset (`StrelkaDemo`) which is used
for testing of workflow.

See also [TinDaisy](https://github.com/ding-lab/tin-daisy) for a workflow wrapper
which includes SomaticWrapper:CWL


## Tags

This was done at some point: 
`docker tag cgc-images.sbgenomics.com/m_wyczalkowski/somatic-wrapper:cwl cgc-images.sbgenomics.com/m_wyczalkowski/somatic-wrapper:cwl-v1.0`
