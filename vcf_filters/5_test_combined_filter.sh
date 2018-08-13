# Testing pyvcf's exensible vcf_flter.py framework

DATAD="/data"  # this works if running inside of docker

STRELKA_VCF="$DATAD/origdata/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="$DATAD/dat/varscan.short.vcf"
VARSCAN_INDEL_VCF="$DATAD/origdata/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/dat/pindel.short.vcf"

VCF=$STRELKA_VCF
OUT="strelka.test.vcf"

CALLER="--caller strelka"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Common configuration file is used for all filters
CONFIG="--config vcf_filter_config.ini"

# Arguments to VAF filter
VAF_FILTER="$MAIN_FILTER --local-script vaf_filter.py"  # filter module
VAF_FILTER_ARGS="vaf --debug $CONFIG $CALLER" 

# Arguments to length filter
LENGTH_FILTER="$MAIN_FILTER --local-script length_filter.py"  # filter module
LENGTH_FILTER_ARGS="indel_length --debug $CONFIG" 

# Arguments to depth filter
DEPTH_FILTER="$MAIN_FILTER --local-script depth_filter.py"  # filter module
DEPTH_FILTER_ARGS="read_depth --debug $CONFIG $CALLER" 

$VAF_FILTER $VCF $VAF_FILTER_ARGS | $LENGTH_FILTER - $LENGTH_FILTER_ARGS | $DEPTH_FILTER - $DEPTH_FILTER_ARGS > $OUT

echo Written to $OUT

