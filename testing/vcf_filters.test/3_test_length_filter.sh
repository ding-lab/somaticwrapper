# Testing length filter using pyvcf's exensible vcf_filter.py framework
source common_config.sh

function run_length_filter {
VCF=$1; shift
XARGS="$@"

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path
FILTER_LOCAL="length_filter.py"  # filter module

$MAIN_FILTER --local-script $FILTER_LOCAL $VCF length $CONFIG $XARGS

}

CONFIG="--config ../../params/vcf_filter_config.ini"

VARSCAN_INDEL_VCF="StrelkaDemo-results/varscan/filter_out/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="StrelkaDemo-results/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf"

run_length_filter $VARSCAN_INDEL_VCF --debug 

#run_length_filter $PINDEL_VCF --bypass


