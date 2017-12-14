#!/bin/bash

# remove DBSNP variants which exist in the COSMIC database
# Creates dbsnp.noCOSMIC.GRCh38.vcf

DATD="/data/image.data/B_Filter"

DBSNP="$DATD/human_9606_b150_GRCh38p7.All.5col.vcf.gz"
COSMICDB="$DATD/CosmicCodingMuts.grch38.v82.vcf.gz"
OUT="$DATD/dbsnp.noCOSMIC.GRCh38.vcf.gz"  

bash ./make_variant_filter.sh $DBSNP $COSMICDB $OUT
