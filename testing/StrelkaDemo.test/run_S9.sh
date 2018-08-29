source project_config.sh

STEP="vep_annotate"

#9 vep_annotate:  
#    --input_vcf s: VCF file to be annotated with vep_annotate.  Required
#    --reference_fasta s: path to reference.  Required
#    --results_dir s: Per-sample analysis results location. Often same as sample name [.] 
#    --assembly s: either "GRCh37" or "GRCh38", used to identify cache file. Optional if not ambigous 
#    --vep_cache_version s: Cache version, e.g. '90', used to identify cache file.  Optional if not ambiguous
#    --vep_cache_dir s: location of VEP cache directory, indicator whether to use online VEP DB lookups.  
#        * if vep_cache_dir is not defined, will perform online VEP DB lookups
#        * If vep_cache_dir is a directory, it indicates location of VEP cache 
#        * If vep_cache_dir is a file ending in .tar.gz, will extract its contents into "./vep-cache" and use VEP cache
#        NOTE: Online VEP database lookups a) uses online database (so cache isn't installed) b) does not use tmp files
#          It is meant to be used for testing and lightweight applications.  Use the cache for better performance.
#          See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
#    --vep_output: Define output format after annotation.  Allowed values: vcf, vep.  [vcf]

INPUT_VCF="results/merged/merged.filtered.vcf"  

ARGS="\
--input_vcf $INPUT_VCF \
--reference_fasta $REFERENCE_FASTA \
--results_dir $RESULTS_DIR \
"

# We rely on online VEP cache lookup for initial testing, so vep_cache_dir is not specified

# --assembly s: either "GRCh37" or "GRCh38", used to identify cache file. Optional if not ambigous 
# --vep_cache_version s: Cache version, e.g. '90', used to identify cache file.  Optional if not ambiguous
# --vep_cache_dir s: location of VEP cache directory, indicator whether to use online VEP DB lookups.  
# --vep_output: Define output format after annotation.  Allowed values: vcf, vep.  [vcf]

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# output: results/vep/output.vcf
