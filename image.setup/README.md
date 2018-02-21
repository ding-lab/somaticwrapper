Scripts here are to be run just once to set up common files (like reference,
etc), which can be used by multiple docker containers.  Scripts must be run from within the docker container.
See [SomaticWrapper.CPTAC3.b1](https://github.com/ding-lab/SomaticWrapper.CPTAC3.b1) for details.

All data is written to the `/image` directory (defined in SomaticWrapper.workflow by
`$IMAGED_H`).  

# Production initialization

Below are the steps to download/generate files needed by SomaticWrapper.  Additional details in
directories `A_Reference`, `B_Filter`, `C_Centromeres`, and `D_VEP`.

## Reference

Download and index a variety of references to `/image/A_Reference`.  Run these scripts as necessary, and more can be readily added.
Alternatively, the reference `.fa` files and associated index files can be simply copied to `/image/A_Reference`.

## Filter

Download and create the variant database. This process consists of three parts:

1. Download dbsnp database
2. Download COSMIC database
3. Remove from dbsnp database those mutations which exist in COSMIC database

This is performed for both the GRCh37 and GRCh38 assemblies.

## Centromeres

Create a centromere BED file.  This is used by Pindel to define regions to be excluded
from analysis for performance reasons.

## VEP

Install the VEP cache for annotating variants.  This is optional, but makes annotation much faster for large datasets.

# Strelka Demo

**TODO**: Update `B_Filter/4_makeStrelkaTestData.sh` so that it writes `dbsnp-StrelkaDemo.noCOSMIC.vcf.gz` to `S_StrelkaDemo/data` and
check it.  Initialization of db will involve copying data to `/image`

Test this on docker
Also, set up docker run script so `/image` and `/import` volumes can be mounted `rw`

Next, describe these steps in the SomaticWrapper.CPTAC3.b1 script


# To run with test Strelka directory

1. Run scripts in `S_StrelkaDemoSetup`
2. 

# To run arbitrary data

1. set up image
2. create config file (example)

