# Set up data directory which will be used by SomaticWrapper
# This script is run outside on the host
# The output is a directory mounted within the container
#
# Copy test data from Strelka for test purposes
# to obtain the strelka data,
#   cd /home/mwyczalk_test/src
#   git clone https://github.com/Illumina/strelka

# Note that we have to either copy BAM/BAI data, or create hard links to it.
# Here we copy. Note that hard links require that the DATA_DIR be in the same 
# partition as the data

# To make VEP annotation work, we need to rename the chrom "demo20" to "20".  This
# is done for both the BAM files and the reference

DATA_DIR="/home/mwyczalk_test/data/docker/data/data"   # this ends up mapping to /data/data within container
REF_DIR="/home/mwyczalk_test/data/docker/data/A_Reference"   # Reference data lives in a different place to keep from confusing SomaticWrapper
STRELKA_DIR="/home/mwyczalk_test/src/strelka"

# Sample Name is used by Somatic Wrapper in the naming of the data files
SAMPLE_NAME="SWtest"

# Create directory structure expected by Somatic Wrapper
SAMPLE_DIR="$DATA_DIR/$SAMPLE_NAME"
mkdir -p $SAMPLE_DIR
mkdir -p $REF_DIR

# define original tumor, normal BAM and index files 
N_BAM="$STRELKA_DIR/src/demo/data/NA12892_demo20.bam"
N_BAI="$STRELKA_DIR/src/demo/data/NA12892_demo20.bam.bai"

T_BAM="$STRELKA_DIR/src/demo/data/NA12891_demo20.bam"
T_BAI="$STRELKA_DIR/src/demo/data/NA12891_demo20.bam.bai"

# Reference data lives in a different place to keep from confusing SomaticWrapper
REF="$STRELKA_DIR/src/demo/data/demo20.fa"
REFI="$STRELKA_DIR/src/demo/data/demo20.fa.fai"


# Normally this is what we'd do:
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
samtools faidx $REF_DIR/demo20.fa

echo writing to $SAMPLE_DIR

# Renaming 'demo20' -> '20' in BAM files and reindexing

samtools view -h $N_BAM | sed 's/demo20/20/g' | samtools view -bT $REF_DIR/demo20.fa - -o $SAMPLE_DIR/${SAMPLE_NAME}.N.bam
samtools index $SAMPLE_DIR/${SAMPLE_NAME}.N.bam $SAMPLE_DIR/${SAMPLE_NAME}.N.bai

samtools view -h $T_BAM | sed 's/demo20/20/g' | samtools view -bT $REF_DIR/demo20.fa - -o $SAMPLE_DIR/${SAMPLE_NAME}.T.bam
samtools index $SAMPLE_DIR/${SAMPLE_NAME}.T.bam $SAMPLE_DIR/${SAMPLE_NAME}.T.bai

echo Written test data to $SAMPLE_DIR and $REF_DIR

