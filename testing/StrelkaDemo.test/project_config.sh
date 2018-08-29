DATAD="/data"
TUMOR_BAM=$DATAD/StrelkaDemoCase.T.bam
NORMAL_BAM=$DATAD/StrelkaDemoCase.N.bam
REFERENCE_FASTA=$DATAD/demo20.fa
STRELKA_CONFIG=$DATAD/strelka.WES.ini
VARSCAN_CONFIG=$DATAD/varscan.WES.ini
PINDEL_CONFIG=$DATAD/pindel.WES.ini
DBSNP_DB=$DATAD/dbsnp-StrelkaDemo.noCOSMIC.vcf.gz
CENTROMERE_BED=$DATAD/ucsc-centromere.GRCh37.bed
#VEP_CACHE_DIR=/home/mwyczalk_test/data/docker/data/D_VEP
ASSEMBLY="GRCh37"
STRELKA_VCF_FILTER_CONFIG="/usr/local/somaticwrapper/vcf_filters/vcf_filter_config.ini"
VARSCAN_VCF_FILTER_CONFIG="/usr/local/somaticwrapper/vcf_filters/vcf_filter_config.ini"
PINDEL_VCF_FILTER_CONFIG="/usr/local/somaticwrapper/vcf_filters/pindel-vcf_filter_config.ini"

RESULTS_DIR="results"
mkdir -p $RESULTS_DIR

SAMPLE="StrelkaDemo"
