# Testing pyvcf's exensible vcf_flter.py framework
# looking at merged data.  Here, pindel samples != (strelka, varscan) samples

DATAD="/data"  # this works if running inside of docker

MERGED_VCF="$DATAD/dat/td2.merged.short.vcf"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
VAF_FILTER_LOCAL="merge_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Arguments to VAF filter
SNV_VAF_ARGS="merged_caller --debug --include strelka-varscan" # --debug"
$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $MERGED_VCF $SNV_VAF_ARGS 

