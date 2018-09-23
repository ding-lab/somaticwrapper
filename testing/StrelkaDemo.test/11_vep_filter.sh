DATAD="/usr/local/somaticwrapper/testing/StrelkaDemo.dat"
source project_config.sh $DATAD

STEP="vep_filter"

# vep_filter:  
#     --input_vcf s: VCF file to be annotated with vep_annotate.  Required
#     --af_filter_config s: configuration file for af (allele frequency) filter.  Required.
#     --classification_filter_config s: configuration file for classification filter.  Required.
#     --bypass_af: Bypass AF filter by retaining all reads
#     --bypass_classification: Bypass Classification filter by retaining all reads
#     --bypass: Bypass all filters
#     --debug: print out processing details to STDERR

INPUT_VCF="results/vep/output_vep.vcf"

# We rely on online VEP cache lookup for StrelkaDemo testing, so vep_cache_dir is not specified
# AF filtering cannot be performed as a result, since that requires cache

# Parameters needed for cache use:
#--vep_cache_dir /image/D_VEP \
#--vep_cache_version 90 \
#--assembly GRCh37 \

ARGS="\
--input_vcf $INPUT_VCF \
--results_dir $RESULTS_DIR \
--af_filter_config $AF_FILTER_CONFIG \
--classification_filter_config $CLASS_FILTER_CONFIG \
--bypass \
"


BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# output: results/vep/output.vcf
