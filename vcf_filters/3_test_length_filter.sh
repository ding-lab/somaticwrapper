# Testing pyvcf's exensible vcf_flter.py framework

DATAD="/data"  # this works if running inside of docker

VARSCAN_VCF="$DATAD/origdata/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/origdata/pindel.out.current_final.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/dat/pindel.short.vcf"

export PYTHONPATH="$DATAD/somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
LENGTH_FILTER_LOCAL="length_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# arguments to length filter
LENGTH_ARGS="indel_length --max_length 10 --debug "

$MAIN_FILTER --local-script $LENGTH_FILTER_LOCAL $VARSCAN_VCF $LENGTH_ARGS


# TESTING
# * pindel - OK
# * varscan indel - OK

