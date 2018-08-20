# Filter merged VCF file to exclude snv calls generated by just one caller.
# Specifically, we exclude calls where "set" INFO field is "strelka" or "varscan".
#   implicitly, pindel, varindel, and strelka-varscan calls are retained
#
# Usage:
#   bash run_merged_filter.sh input.vcf output.vcf [args ...]
# args is optional argument passed to all filters, e.g., --debug

VCF=$1; shift
OUT=$2; shift
XARG=$@  # https://stackoverflow.com/questions/1537673/how-do-i-forward-parameters-to-other-command-in-bash-script

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"

MERGE_FILTER="vcf_filter.py --no-filtered --local-script merge_filter.py"  # filter module
MERGE_FILTER_ARGS="merged_caller --exclude strelka,varscan $XARG " 

$MERGE_FILTER $VCF $MERGE_FILTER_ARGS > $OUT

>&2 echo Written to $OUT

