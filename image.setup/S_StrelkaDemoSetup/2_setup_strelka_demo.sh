# Set up data directory which will be used by SomaticWrapper
#
# Using Strelka test data from https://github.com/Illumina/strelka/tree/master/src/demo/data
# This is a compact dataset which has variants on chromosome "demo20"
# Using NA12892_demo20.bam as Normal
#       NA12891_demo20.bam as Tumor

# To make VEP annotation work, we need to rename the chrom "demo20" to "20".  This
# is done for both the BAM files and the reference. 

# SomaticWrapper directory structure
# /image/A_Reference - reference data
# /image/A_Reference/demo20.fa - test reference
# /data/SWtest - per-run SomaticWrapper results
# /import/StrelkaTestData/SWtest.N.bam, SWtest.T.bam - tumor and normal BAMs

DATAD_C="/import/StrelkaTestData" # Where test BAMs live
REFDIR_C="/image/A_Reference"       # Reference data lives in a different place to keep from confusing SomaticWrapper

# Sample Name (aka Run Name) is used by Somatic Wrapper in the naming of the data files
# TODO: This usage of "sample name" is inconsistent with that used by SomaticWrapper.Workflow, where SampleName refers
# to a single BAM file.  This should be fixed or at least noted
SAMPLE_NAME="SWtest"

# Create directory structure expected by Somatic Wrapper.  
mkdir -p $DATAD_C
mkdir -p $REFDIR_C

# define original tumor, normal BAM and index files 
N_BAM="$IMPORTD_C/NA12892_demo20.bam"
N_BAI="$IMPORTD_C/NA12892_demo20.bam.bai"

T_BAM="$IMPORTD_C/NA12891_demo20.bam"
T_BAI="$IMPORTD_C/NA12891_demo20.bam.bai"

# Reference data lives in a different place to keep from confusing SomaticWrapper
REF="$IMAGED_C/demo20.fa"
REFI="$IMAGED_C/demo20.fa.fai"

#####
# so that this example works with vep annotation (which assumes real chrom names), we need to rename chrom "demo20" to standard name "20"
# fix the reference first
sed 's/demo20/20/g' $REF > $REFDIR_C/demo20.fa
# and index it
bash ../A_Reference/prepare_reference.sh $REFDIR_C/demo20.fa

# Renaming 'demo20' -> '20' in BAM files and reindexing
# 
samtools view -h $N_BAM | sed 's/demo20/20/g' | samtools view -bT $REFDIR_C/demo20.fa - -o $DATAD_C/${SAMPLE_NAME}.N.bam
samtools index $DATAD_C/${SAMPLE_NAME}.N.bam $DATAD_C/${SAMPLE_NAME}.N.bai

samtools view -h $T_BAM | sed 's/demo20/20/g' | samtools view -bT $REFDIR_C/demo20.fa - -o $DATAD_C/${SAMPLE_NAME}.T.bam
samtools index $DATAD_C/${SAMPLE_NAME}.T.bam $DATAD_C/${SAMPLE_NAME}.T.bai

echo Written test data to $DATAD_C and $REFDIR_C

# Copy demo dbsnp to /image directory
# This file is created by ../B_Filter/4_makeStrelkaTestData.sh
# Here, we just copy it to the appropriate directory
mkdir -p /image/B_Filter
cp data/dbsnp-StrelkaDemo.noCOSMIC.vcf.gz /image/B_Filter
