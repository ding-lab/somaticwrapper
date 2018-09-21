# run VAF, length, and read depth filters on a VCF
# Usage:
#   bash run_combined_filter.sh input.vcf caller config.ini output.vcf [args]
# where CALLER is one of strelka, varscan, or pindel
# config.ini is configuration file used by all filters
# args are zero or more optional arguments:
#   --bypass_vaf, --bypass_depth, --bypass_length will skip just that filter
#   --bypass will skip all filters
#   --debug will print out debug info to STDERR for all filters
#   --debug_vaf, debug_depth, debug_length - debug specific to filter
# If output.vcf is -, write to stdout

if [ "$#" -lt 4 ]; then
    >&2 echo Error: Wrong number of arguments
    >&2 echo Usage:  bash run_combined_filter.sh input.vcf caller config.ini output.vcf [args]

    exit 1
fi


VCF=$1 ; shift
CALLER=$1 ; shift
CONFIG_FN=$1 ; shift
OUT=$1 ; shift
XARG="$@"  # optional arguments

CALLER_ARG="--caller $CALLER"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"

# parse XARG to catch various bypass and debug options.
# --bypass will bypass all
# if this is not a bypass or debug arg then add it to both filters
for ARG in $XARG; do
    if [ "$ARG" == "--bypass_vaf" ]; then
        VAF_ARG="$VAF_ARG --bypass"
    elif [ "$ARG" == "--bypass_length" ]; then
        LENGTH_ARG="$LENGTH_ARG --bypass"
    elif [ "$ARG" == "--bypass_depth" ]; then
        DEPTH_ARG="$DEPTH_ARG --bypass"
    elif [ "$ARG" == "--bypass" ]; then
        VAF_ARG="$VAF_ARG --bypass"
        LENGTH_ARG="$LENGTH_ARG --bypass"
        DEPTH_ARG="$DEPTH_ARG --bypass"
    elif [ "$ARG" == "--debug_vaf" ]; then
        VAF_ARG="$VAF_ARG --debug"
    elif [ "$ARG" == "--debug_length" ]; then
        LENGTH_ARG="$LENGTH_ARG --debug"
    elif [ "$ARG" == "--debug_depth" ]; then
        DEPTH_ARG="$DEPTH_ARG --debug"
    elif [ "$ARG" == "--debug" ]; then
        VAF_ARG="$VAF_ARG --debug"
        LENGTH_ARG="$LENGTH_ARG --debug"
        DEPTH_ARG="$DEPTH_ARG --debug"
    else
        VAF_ARG="$VAF_ARG $ARG"
        LENGTH_ARG="$LENGTH_ARG $ARG"
        DEPTH_ARG="$DEPTH_ARG $ARG"
    fi
done


MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Common configuration file is used for all filters
CONFIG="--config $CONFIG_FN"

# Arguments to VAF filter
VAF_FILTER="vcf_filter.py --no-filtered --local-script vaf_filter.py"  # filter module
VAF_FILTER_ARGS="vaf $VAF_ARG $CONFIG $CALLER_ARG" 

# Arguments to length filter
LENGTH_FILTER="vcf_filter.py --no-filtered --local-script length_filter.py"  # filter module
LENGTH_FILTER_ARGS="length $LENGTH_ARG $CONFIG" 

# Arguments to depth filter
DEPTH_FILTER="vcf_filter.py --no-filtered --local-script depth_filter.py"  # filter module
DEPTH_FILTER_ARGS="read_depth $DEPTH_ARG $CONFIG $CALLER_ARG" 

if [ $OUT == '-' ]; then

$VAF_FILTER $VCF $VAF_FILTER_ARGS | $LENGTH_FILTER - $LENGTH_FILTER_ARGS | $DEPTH_FILTER - $DEPTH_FILTER_ARGS 

else

$VAF_FILTER $VCF $VAF_FILTER_ARGS | $LENGTH_FILTER - $LENGTH_FILTER_ARGS | $DEPTH_FILTER - $DEPTH_FILTER_ARGS > $OUT

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

