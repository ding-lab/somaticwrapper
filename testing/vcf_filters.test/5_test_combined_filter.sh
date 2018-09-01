# Testing combined filter (
source common_config.sh

function run_combined_filter {
CALLER=$1; shift
VCF=$1; shift
OUT=$1; shift
XARGS="$@"

RUN="../../src/vcf_filters/run_combined_vcf_filter.sh"
bash $RUN $VCF $CALLER $CONFIG $OUT $XARGS

}

CONFIG="../../params/vcf_filter_config.ini"

STRELKA_VCF="$DATAD/root/s3_parse_strelka/results/strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="$DATAD/root/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
VARSCAN_INDEL_VCF="$DATAD/root/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/root/s7_parse_pindel/results/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf"


run_combined_filter strelka $STRELKA_VCF strelka.tmp.vcf 

run_combined_filter varscan $VARSCAN_VCF varscan.tmp.vcf --debug

run_combined_filter varscan $VARSCAN_INDEL_VCF varindel.tmp.vcf --bypass --debug

run_combined_filter pindel $PINDEL_VCF pindel.tmp.vcf --bypass 

# ORIGINAL below
exit

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



