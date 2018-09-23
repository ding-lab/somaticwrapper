# run Allele Frequency and Consequence filters on a VCF
# Usage:
#   bash run_combined_af_classification_filter.sh input.vcf af_config.ini classification_config.ini output.vcf [args]
# args are zero or more optional arguments:
#   --bypass_af, --bypass_classification will skip just that filter
#   --bypass will skip all filters
#   --debug will print out debug info to STDERR
#   --debug_af and --debug_classification specific to that filter

VCF=$1 ; shift
AF_CONFIG=$1 ; shift
CLASS_CONFIG=$1 ; shift
OUT=$1 ; shift
XARG="$@"  # optional argument passed to all filters, may be --debug

# parse XARG to catch various bypass and debug options.
# --bypass will bypass both
# --debug will debug both
# if this is not a bypass arg then add it to both filters 
for ARG in $XARG; do
    if [ "$ARG" == "--bypass_af" ]; then
        AF_ARG="$AF_ARG --bypass"
    elif [ "$ARG" == "--bypass_classification" ]; then
        CLASS_ARG="$CLASS_ARG --bypass"
    elif [ "$ARG" == "--bypass" ]; then
        AF_ARG="$AF_ARG --bypass"
        CLASS_ARG="$CLASS_ARG --bypass"
    elif [ "$ARG" == "--debug_af" ]; then
        AF_ARG="$AF_ARG --debug"
    elif [ "$ARG" == "--debug_classification" ]; then
        CLASS_ARG="$CLASS_ARG --debug"
    elif [ "$ARG" == "--debug" ]; then
        AF_ARG="$AF_ARG --debug"
        CLASS_ARG="$CLASS_ARG --debug"
    else
        CLASS_ARG="$CLASS_ARG $ARG"
        AF_ARG="$AF_ARG $ARG"
    fi
done

export PYTHONPATH="/usr/local/somaticwrapper/src/vcf_filters:$PYTHONPATH"

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Arguments to AF filter
# --bypass_if_missing will keep filter from exiting if AF_MAX is missing from e.g. VEP DB annotation
AF_FILTER="vcf_filter.py --no-filtered --local-script af_filter.py"  # filter module
AF_FILTER_ARGS="af $AF_ARG --config $AF_CONFIG --input_vcf $VCF --bypass_if_missing"

# Arguments to classification filter
CLASS_FILTER="vcf_filter.py --no-filtered --local-script classification_filter.py"  # filter module
CLASS_FILTER_ARGS="classification $CLASS_ARG --config $CLASS_CONFIG --input_vcf $VCF" 

if [ $OUT == '-' ]; then

$AF_FILTER $VCF $AF_FILTER_ARGS | $CLASS_FILTER - $CLASS_FILTER_ARGS 

else

$AF_FILTER $VCF $AF_FILTER_ARGS | $CLASS_FILTER - $CLASS_FILTER_ARGS > $OUT

fi

# Evaluate return value for chain of pipes; see https://stackoverflow.com/questions/90418/exit-shell-script-based-on-process-exit-code
rcs=${PIPESTATUS[*]};
for rc in ${rcs}; do
    if [[ $rc != 0 ]]; then
        >&2 echo Fatal error.  Exiting.
        exit $rc;
    fi;
done

>&2 echo Written to $OUT

