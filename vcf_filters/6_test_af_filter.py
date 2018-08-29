# Testing depth filter using pyvcf's exensible vcf_filter.py framework

DATAD="/data"  # this works if running inside of docker

VEP_VCF="../testing/C3N-01649.test/C3N-01649.test.vcf"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
AF_FILTER_LOCAL="af_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# arguments to depth filter
ARGS="af --debug"

$MAIN_FILTER --local-script $AF_FILTER_LOCAL $VEP_VCF $ARGS 
