export TUMOR_BAM=/data/StrelkaDemoCase.T.bam
export NORMAL_BAM=/data/StrelkaDemoCase.N.bam
export REFERENCE_FASTA=/data/demo20.fa
export STRELKA_CONFIG=/data/strelka.WES.ini
export DBSNP_DB=/data/dbsnp-StrelkaDemo.noCOSMIC.vcf.gz
export CENTROMERE_BED=/data/ucsc-centromere.GRCh37.bed
#export VEP_CACHE_DIR=/home/mwyczalk_test/data/docker/data/D_VEP

export OUTPUT_DIR=./results  #/diskmnt/Projects/TinDaisy/results

export VARSCAN_CONFIG=../../params/varscan.WES.ini

# try to have all output go to output_dir
mkdir -p $OUTPUT_DIR

ARGS="--tumor_bam $TUMOR_BAM \
--normal_bam $NORMAL_BAM \
--reference_fasta $REFERENCE_FASTA \
--strelka_config $STRELKA_CONFIG \
--is_strelka2 1"

perl ../../SomaticWrapper.pl $ARGS --results_dir $OUTPUT_DIR 1
