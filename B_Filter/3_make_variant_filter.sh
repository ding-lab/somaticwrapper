#!/bin/bash

# remove DBSNP variants which exist in the COSMIC database
# Creates dbsnp.noCOSMIC.vcf

DATD="/data/B_Filter"
JAR="/usr/local/snpEff/SnpSift.jar"

#DBSNP="$DATD/short.5col.vcf.gz"
DBSNP="$DATD/human_9606_b150_GRCh37p13.5col.vcf.gz"
COSMICDB="$DATD/CosmicCodingMuts.grch37.v82.vcf.gz"

# Output VCF
OUT="$DATD/dbsnp.noCOSMIC.vcf"  # equivalant to pass.cosmic.vcf

# if ANNO is defined, write out this intermediate file.  useful for debugging
# if ANNO is undefined (by commenting it out), will use pipeline vesion which does not generate an
#   intermediate file.  

#ANNO="$DATD/dbsnp_cosmic_anno.vcf"

if [ ! -z "$ANNO" ]; then

echo Writing intermediate $ANNO

java -jar $JAR annotate -id $COSMICDB $DBSNP  > $ANNO
echo Written to $ANNO

java -jar $JAR filter -n " (exists ID) & (ID =~ 'COSM' ) " -f $ANNO > $OUT
echo Written to $OUT

else

java -jar $JAR annotate -id $COSMICDB $DBSNP  | java -jar $JAR filter -n " (exists ID) & (ID =~ 'COSM' ) "  > $OUT
echo Written to $OUT

fi


