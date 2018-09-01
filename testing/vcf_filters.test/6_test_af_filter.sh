# Testing AF filter using pyvcf's exensible vcf_filter.py framework
source common_config.sh

# Note: require the --flag_pick flag when running vep
VEP_VCF="../C3N-01649.test/C3N-01649.results/vep/output.vcf"
CONFIG="--config ../../params/af_filter_config.ini"

AF_FILTER_LOCAL="af_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# arguments to depth filter
ARGS="af --debug --input_vcf $VEP_VCF --bypass"

$MAIN_FILTER --local-script $AF_FILTER_LOCAL $VEP_VCF $ARGS $CONFIG
