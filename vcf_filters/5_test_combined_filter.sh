# Testing pyvcf's exensible vcf_flter.py framework

DATAD="/data"  # this works if running inside of docker

STRELKA_VCF="$DATAD/origdata/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="$DATAD/dat/varscan.short.vcf"
VARSCAN_INDEL_VCF="$DATAD/origdata/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/dat/pindel.short.vcf"

VCF=$STRELKA_VCF
CALLER="strelka"
CONFIG="vcf_filter_config.ini"
OUT="strelka.test.vcf"

bash run_combined_filter.sh $VCF $CALLER $CONFIG $OUT 
