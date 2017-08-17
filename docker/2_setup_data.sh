# Set up data directory which will be used by SomaticWrapper
# This will be a directory mounted within the container

# Note that we have to either copy BAM/BAI data, or create hard links to it.
# Here we create hard links, but note that this requires that the DATA_DIR be in the same 
# partition as the data

DATA_DIR="/home/mwyczalk_test/data/docker/data/data"

# Sample Name is used by Somatic Wrapper in the naming of the data files
SAMPLE_NAME="SWtest"

# define tumor, normal BAM and index files 
N_BAM="/home/mwyczalk_test/src/SomaticWrapper/strelka/src/demo/data/NA12892_dupmark_chr20_region.bam"
N_BAI="/home/mwyczalk_test/src/SomaticWrapper/strelka/src/demo/data/NA12892_dupmark_chr20_region.bam.bai"

T_BAM="/home/mwyczalk_test/src/SomaticWrapper/strelka/src/demo/data/NA12891_dupmark_chr20_region.bam"
T_BAI="/home/mwyczalk_test/src/SomaticWrapper/strelka/src/demo/data/NA12891_dupmark_chr20_region.bam.bai"

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



