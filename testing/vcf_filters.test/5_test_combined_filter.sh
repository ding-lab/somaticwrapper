# Testing pyvcf's exensible vcf_filter.py framework

DATAD="/data"  # this works if running inside of docker

STRELKA_VCF="$DATAD/origdata/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="$DATAD/dat/varscan.short.vcf"
VARSCAN_INDEL_VCF="$DATAD/origdata/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/dat/pindel.short.vcf"

CONFIG="../params/vcf_filter_config.ini"

VCF=$STRELKA_VCF
CALLER="strelka"
OUT="strelka.test.vcf"
bash run_combined_filter.sh $VCF $CALLER $CONFIG $OUT 

VCF=$VARSCAN_VCF
CALLER="varscan"
OUT="varscan.test.vcf"
bash run_combined_filter.sh $VCF $CALLER $CONFIG $OUT 

VCF=$VARSCAN_INDEL_VCF
CALLER="varscan"
OUT="varindel.test.vcf"
bash run_combined_filter.sh $VCF $CALLER $CONFIG $OUT 

VCF=$PINDEL_VCF
CALLER="pindel"
OUT="pindel.test.vcf"
CONFIG="../params/pindel-vcf_filter_config.ini"
bash run_combined_filter.sh $VCF $CALLER $CONFIG $OUT 



