# pass root of data directory as first argument, or use /data by default
if [ -z $1 ]; then
    DATAD="/data"
else
    DATAD=$1
fi
PARAMS="/usr/local/somaticwrapper/params"

TUMOR_BAM="$DATAD/StrelkaDemoCase.T.bam"
NORMAL_BAM="$DATAD/StrelkaDemoCase.N.bam"
REFERENCE_FASTA="$DATAD/demo20.fa"
STRELKA_CONFIG="$PARAMS/strelka.WES.ini"
VARSCAN_CONFIG="$PARAMS/varscan.WES.ini"
PINDEL_CONFIG="$PARAMS/pindel.WES.ini"
DBSNP_DB="$DATAD/dbsnp-StrelkaDemo.noCOSMIC.vcf.gz"
CENTROMERE_BED="$DATAD/ucsc-centromere.GRCh37.bed"
#VEP_CACHE_DIR=/home/mwyczalk_test/data/docker/data/D_VEP
ASSEMBLY="GRCh37"

STRELKA_VCF_FILTER_CONFIG="$PARAMS/vcf_filter_config-strelka.ini"
VARSCAN_VCF_FILTER_CONFIG="$PARAMS/vcf_filter_config-varscan.ini"
PINDEL_VCF_FILTER_CONFIG="$PARAMS/vcf_filter_config-pindel.ini"
AF_FILTER_CONFIG="$PARAMS/af_filter_config.ini"
CLASS_FILTER_CONFIG="$PARAMS/classification_filter_config.ini"

RESULTS_DIR="results"
mkdir -p $RESULTS_DIR

SAMPLE="StrelkaDemo"
