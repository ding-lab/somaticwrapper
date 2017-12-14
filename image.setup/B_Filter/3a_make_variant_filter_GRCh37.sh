#!/bin/bash

# remove DBSNP variants which exist in the COSMIC database
# Creates dbsnp.noCOSMIC.GRCh37.vcf

DATD="/data/image.data/B_Filter"

DBSNP="$DATD/human_9606_b150_GRCh37p13.All.5col.vcf.gz"
COSMICDB="$DATD/CosmicCodingMuts.grch37.v82.vcf.gz"
OUT="$DATD/dbsnp.noCOSMIC.GRCh37.vcf.gz"  

bash ./make_variant_filter.sh $DBSNP $COSMICDB $OUT
