# Set up data directory which will be used by SomaticWrapper
#
# Using Strelka test data from https://github.com/Illumina/strelka/tree/master/src/demo/data
# This is a compact dataset which has variants on chromosome "demo20"
# Using NA12892_demo20.bam as Normal
#       NA12891_demo20.bam as Tumor

# To make VEP annotation work, we need to rename the chrom "demo20" to "20".  This
# is done for both the BAM files and the reference. 

# SomaticWrapper directory structure
# /data/A_Reference - reference data
# /data/A_Reference/demo20.fa - test reference
# /data/data/SWtest - corresponds to individual sample SWtest for analysis (SAMPLE_DIR) 
# /data/data/SWtest.N.bam, SWtest.T.bam - tumor and normal BAMs (may be links)

DATA_BASE="/data"   
DATD="$DATA_BASE/S_StrelkaTestData" # Where test data lives (BAM and reference)
REF_DIR="$DATA_BASE/A_Reference"       # Reference data lives in a different place to keep from confusing SomaticWrapper

STRELKA_DIR="/home/mwyczalk_test/src/strelka"

# Sample Name is used by Somatic Wrapper in the naming of the data files
SAMPLE_NAME="SWtest"

# Create directory structure expected by Somatic Wrapper
SAMPLE_DIR="$DATA_BASE/data/$SAMPLE_NAME"
mkdir -p $SAMPLE_DIR
mkdir -p $REF_DIR

# define original tumor, normal BAM and index files 
N_BAM="$DATD/NA12892_demo20.bam"
N_BAI="$DATD/NA12892_demo20.bam.bai"

T_BAM="$DATD/NA12891_demo20.bam"
T_BAI="$DATD/NA12891_demo20.bam.bai"

# Reference data lives in a different place to keep from confusing SomaticWrapper
REF="$DATD/demo20.fa"
REFI="$DATD/demo20.fa.fai"


# If not remapping chrom, we would simply copy or make links to original data as,
#cp $REF $REF_DIR
#cp $REFI $REF_DIR

#cp $N_BAM $SAMPLE_DIR/${SAMPLE_NAME}.N.bam
#cp $N_BAI $SAMPLE_DIR/${SAMPLE_NAME}.N.bai
#cp $T_BAM $SAMPLE_DIR/${SAMPLE_NAME}.T.bam
#cp $T_BAI $SAMPLE_DIR/${SAMPLE_NAME}.T.bai

#####
# so that this example works with vep annotation (which assumes real chrom names), we need to rename chrom "demo20" to standard name "20"
# fix the reference first
sed 's/demo20/20/g' $REF > $REF_DIR/demo20.fa
# and index it
bash ../A_Reference/prepare_reference.sh $REF_DIR/demo20.fa

# Renaming 'demo20' -> '20' in BAM files and reindexing
# 
samtools view -h $N_BAM | sed 's/demo20/20/g' | samtools view -bT $REF_DIR/demo20.fa - -o $SAMPLE_DIR/${SAMPLE_NAME}.N.bam
samtools index $SAMPLE_DIR/${SAMPLE_NAME}.N.bam $SAMPLE_DIR/${SAMPLE_NAME}.N.bai

samtools view -h $T_BAM | sed 's/demo20/20/g' | samtools view -bT $REF_DIR/demo20.fa - -o $SAMPLE_DIR/${SAMPLE_NAME}.T.bam
samtools index $SAMPLE_DIR/${SAMPLE_NAME}.T.bam $SAMPLE_DIR/${SAMPLE_NAME}.T.bai

echo Written test data to $SAMPLE_DIR and $REF_DIR

