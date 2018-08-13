
# run VAF, length, and read depth filters on a VCF
# Usage:
#   bash run_combined_filter.sh input.vcf CALLER config.ini output.vcf [args]
# where CALLER is one of strelka, varscan, or pindel
# config.ini is configuration file used by all filters
# args is optional argument passed to all filters, e.g., --debug

VCF=$1
CALLER=$2
CONFIG_FN=$3
OUT=$4
XARG=$5  # optional argument passed to all filters, may be --debug

CALLER_ARG="--caller $CALLER"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Common configuration file is used for all filters
CONFIG="--config $CONFIG_FN"

# Arguments to VAF filter
VAF_FILTER="vcf_filter.py --no-filtered --local-script vaf_filter.py"  # filter module
VAF_FILTER_ARGS="vaf $XARG $CONFIG $CALLER_ARG" 

# Arguments to length filter
LENGTH_FILTER="vcf_filter.py --no-filtered --local-script length_filter.py"  # filter module
LENGTH_FILTER_ARGS="indel_length $XARG $CONFIG" 

# Arguments to depth filter
DEPTH_FILTER="vcf_filter.py --no-filtered --local-script depth_filter.py"  # filter module
DEPTH_FILTER_ARGS="read_depth $XARG $CONFIG $CALLER_ARG" 

$VAF_FILTER $VCF $VAF_FILTER_ARGS | $LENGTH_FILTER - $LENGTH_FILTER_ARGS | $DEPTH_FILTER - $DEPTH_FILTER_ARGS > $OUT

>&2 echo Written to $OUT

