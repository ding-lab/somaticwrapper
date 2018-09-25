# Testing depth filter using pyvcf's exensible vcf_filter.py framework
source common_config.sh

function run_depth_filter {
VCF=$1; shift
CONFIG=$1; shift
XARGS="$@"

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path
FILTER_LOCAL="depth_filter.py"  # filter module

$MAIN_FILTER --local-script $FILTER_LOCAL $VCF read_depth $CONFIG $XARGS

}

STRELKA_CONFIG="--config ../../params/vcf_filter_config-strelka.ini"
VARSCAN_CONFIG="--config ../../params/vcf_filter_config-varscan.ini"
PINDEL_CONFIG="--config ../../params/vcf_filter_config-pindel.ini"

STRELKA_VCF="StrelkaDemo-results/strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="StrelkaDemo-results/varscan/filter_out/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
VARSCAN_INDEL_VCF="StrelkaDemo-results/varscan/filter_out/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="StrelkaDemo-results/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf"


run_depth_filter $STRELKA_VCF $STRELKA_CONFIG --debug

run_depth_filter $VARSCAN_VCF $VARSCAN_CONFIG --bypass --debug

run_depth_filter $VARSCAN_INDEL_VCF $VARSCAN_CONFIG --bypass --debug

run_depth_filter $PINDEL_VCF $PINDEL_CONFIG --bypass --debug
