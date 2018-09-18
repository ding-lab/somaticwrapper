# Testing combined filter (
source common_config.sh

function run_combined_filter {
CALLER=$1; shift
VCF=$1; shift
OUT=$1; shift
XARGS="$@"

RUN="../../src/vcf_filters/run_vaf_length_depth_filters.sh"
bash $RUN $VCF $CALLER $CONFIG $OUT $XARGS

}

CONFIG="../../params/vcf_filter_config.ini"

STRELKA_VCF="StrelkaDemo-results/strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="StrelkaDemo-results/varscan/filter_out/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
VARSCAN_INDEL_VCF="StrelkaDemo-results/varscan/filter_out/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="StrelkaDemo-results/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf"

mkdir -p tmp

run_combined_filter strelka $STRELKA_VCF - --debug

run_combined_filter varscan $VARSCAN_VCF tmp/varscan.tmp.vcf --bypass_vaf
run_combined_filter varscan $VARSCAN_VCF tmp/varscan2.tmp.vcf --bypass_depth
run_combined_filter varscan $VARSCAN_VCF tmp/varscan3.tmp.vcf --bypass_depth --bypass

run_combined_filter varscan $VARSCAN_INDEL_VCF tmp/varindel.tmp.vcf --bypass

run_combined_filter pindel $PINDEL_VCF tmp/pindel.tmp.vcf --bypass_length


