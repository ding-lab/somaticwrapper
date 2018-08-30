
# run Allele Frequency and Consequence filters on a VCF
# Usage:
#   bash run_combined_af_classification_filter.sh input.vcf af_config.ini classification_config.ini output.vcf [args]
# args is optional argument passed to all filters, e.g., --debug

VCF=$1
AF_CONFIG=$2
CLASS_CONFIG=$3
OUT=$4
XARG=$5  # optional argument passed to all filters, may be --debug

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Arguments to AF filter
AF_FILTER="vcf_filter.py --no-filtered --local-script af_filter.py"  # filter module
AF_FILTER_ARGS="af $XARG --config $AF_CONFIG --input_vcf $VCF" 

# Arguments to classification filter
CLASS_FILTER="vcf_filter.py --no-filtered --local-script classification_filter.py"  # filter module
CLASS_FILTER_ARGS="classification $XARG --config $CLASS_CONFIG --input_vcf $VCF" 


$AF_FILTER $VCF $AF_FILTER_ARGS | $CLASS_FILTER - $CLASS_FILTER_ARGS > $OUT

>&2 echo Written to $OUT

