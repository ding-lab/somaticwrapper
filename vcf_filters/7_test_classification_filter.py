# Testing depth filter using pyvcf's exensible vcf_filter.py framework

DATAD="/data"  # this works if running inside of docker

# Note: require the --flag_pick flag when running vep
VEP_VCF="../testing/C3N-01649.test/C3N-01649.results/vep/output.vcf"
CONFIG="--config classification_filter_config.ini"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
AF_FILTER_LOCAL="classification_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# arguments to depth filter
ARGS="classification --input_vcf $VEP_VCF "

$MAIN_FILTER --local-script $AF_FILTER_LOCAL $VEP_VCF $ARGS $CONFIG
