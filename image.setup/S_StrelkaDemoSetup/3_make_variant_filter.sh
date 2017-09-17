#!/bin/bash

# remove DBSNP variants which exist in the COSMIC database
# Creates dbsnp.noCOSMIC.vcf

DATD="/data/B_Filter"
mkdir -p $DATD
JAR="/usr/local/snpEff/SnpSift.jar"

# Using VCFs distributed with this project.  These were created by B_Filter/4_makeStrelkaTestData.sh
# dbSnP from human_9606_b150_GRCh37p13 00-All.vcf
# Cosmic from grch37 version v82
DBSNP="data/dbSnP.vcf.gz"
COSMICDB="data/COSMIC.vcf.gz"

# if ANNO is defined, write out this intermediate file and don't compress final VCF.  useful for debugging
# if ANNO is undefined (by commenting it out), will use pipeline vesion which does not generate an
#   intermediate file.  also compresses and indexes resulting file

#ANNO="$DATD/dbsnp_cosmic_anno.vcf"

# Output VCF.  Append '.gz' if ANNO is not defined
OUT="$DATD/dbsnp.noCOSMIC.vcf.gz"  

if [ ! -z "$ANNO" ]; then

echo Writing intermediate $ANNO

java -jar $JAR annotate -id $COSMICDB $DBSNP  > $ANNO
echo Written to $ANNO

java -jar $JAR filter -n " (exists ID) & (ID =~ 'COSM' ) " -f $ANNO > $OUT
echo Written to $OUT

else

java -jar $JAR annotate -id $COSMICDB $DBSNP  | java -jar $JAR filter -n " (exists ID) & (ID =~ 'COSM' ) " | bgzip > $OUT

tabix -p vcf $OUT
echo Written and indexed $OUT

fi


