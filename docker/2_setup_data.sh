# Set up data directory which will be used by SomaticWrapper
# This will be a directory mounted within the container
# Copy test data from Strelka for test purposes
# to obtain the strelka data,
#   cd /home/mwyczalk_test/src
#   git clone https://github.com/Illumina/strelka

# Note that we have to either copy BAM/BAI data, or create hard links to it.
# Here we copy. Note that hard links require that the DATA_DIR be in the same 
# partition as the data

DATA_DIR="/home/mwyczalk_test/data/docker/data/data"
STRELKA_DIR="/home/mwyczalk_test/src/strelka"

# Sample Name is used by Somatic Wrapper in the naming of the data files
SAMPLE_NAME="SWtest"

# define tumor, normal BAM and index files 
N_BAM="$STRELKA_DIR/src/demo/data/NA12892_demo20.bam"
N_BAI="$STRELKA_DIR/src/demo/data/NA12892_demo20.bam.bai"

T_BAM="$STRELKA_DIR/src/demo/data/NA12891_demo20.bam"
T_BAI="$STRELKA_DIR/src/demo/data/NA12891_demo20.bam.bai"

# also copy reference
REF="$STRELKA_DIR/src/demo/data/demo20.fa"
REFI="$STRELKA_DIR/src/demo/data/demo20.fa.fai"

# Create directory structure expected by Somatic Wrapper
mkdir -p $DATA_DIR/$SAMPLE_NAME

cd $DATA_DIR/$SAMPLE_NAME

# Make hard links to data
#ln -i $N_BAM ${SAMPLE_NAME}.N.bam
#ln -i $N_BAI ${SAMPLE_NAME}.N.bai
#ln -i $T_BAM ${SAMPLE_NAME}.T.bam
#ln -i $T_BAI ${SAMPLE_NAME}.T.bai

# Copy data
cp $N_BAM ${SAMPLE_NAME}.N.bam
cp $N_BAI ${SAMPLE_NAME}.N.bai
cp $T_BAM ${SAMPLE_NAME}.T.bam
cp $T_BAI ${SAMPLE_NAME}.T.bai

# Reference data lives in a different place to keep from confusing SomaticWrapper
mkdir -p $DATA_DIR/../A_Reference
cp $REF $DATA_DIR/../A_Reference
cp $REFI $DATA_DIR/../A_Reference

echo Written test data to $DATA_DIR/$SAMPLE_NAME and $DATA_DIR/../A_Reference

