DATAD="/usr/local/somaticwrapper/testing/StrelkaDemo.dat"
source project_config.sh $DATAD

STEP="dbsnp_filter"

# dbsnp_filter:
#     --dbsnp_db s: database for dbSNP filtering.  Step will be skipped if not defined
#     --input_vcf s: VCF file to process.  Required
#     --bypass_dbsnp: Apply dbSnP annotation to VCF with no further filtering
#     --bypass: Same as --bypass_dbsnp
#     --debug: print out processing details to STDERR

INPUT_VCF="results/merged/merged.filtered.vcf"

# We rely on online VEP cache lookup for StrelkaDemo testing, so vep_cache_dir is not specified
# AF filtering cannot be performed as a result, since that requires cache

# Parameters needed for cache use:
#--vep_cache_dir /image/D_VEP \
#--vep_cache_version 90 \
#--assembly GRCh37 \

ARGS="\
--input_vcf $INPUT_VCF \
--results_dir $RESULTS_DIR \
--dbsnp_db $DBSNP_DB \
"
#--bypass \


BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# output: results/vep/output.vcf
