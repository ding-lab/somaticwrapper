#!/bin/bash

# remove DBSNP variants which exist in the COSMIC database
# Creates dbsnp.noCOSMIC.vcf

JAR="/usr/local/snpEff/SnpSift.jar"

# Usage: make_variant_filter.sh DBSNP COSMICDB OUT
# example:
# DBSNP="$DATD/human_9606_b150_GRCh37p13.All.5col.vcf.gz"
# COSMICDB="$DATD/CosmicCodingMuts.grch37.v82.vcf.gz"
# OUT="$DATD/dbsnp.noCOSMIC.GRCh37.vcf.gz"  

DBSNP=$1
COSMICDB=$2
OUT=$3

if [ -e $OUT ]; then
	echo $OUT exists.  Please delete before running this again.
	exit
fi

echo Processing $DBSNP and $COSMICDB
echo Writing to $OUT

# if ANNO is defined, write out this intermediate file and don't compress final VCF.  useful for debugging
# if ANNO is undefined (by commenting it out), will use pipeline vesion which does not generate an
#   intermediate file.  also compresses and indexes resulting file

#ANNO="$DATD/dbsnp_cosmic_anno.vcf"

if [ ! -z "$ANNO" ]; then

# Remove .gz extension if it exists, since not saving compressed file here
OUT=${OUT%.gz}

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


