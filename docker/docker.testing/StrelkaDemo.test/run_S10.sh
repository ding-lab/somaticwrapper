source project_config.sh

# This step is not implemented for StrelkaDemo because 
# vcf2maf.pl requires a VEP cache, which is in general not installed for the demo
# With a cache available (either installed or as .tar.gz) 

STEP="vcf_2_maf"
#10 vcf_2_maf:
#    --input_vcf s: VCF file to be annotated with vep_annotate.  Required
#    --reference_fasta s: path to reference.  Required
#    --results_dir s: Per-sample analysis results location. Often same as sample name [.] 
#    --assembly s: either "GRCh37" or "GRCh38", used to identify cache file. Optional if not ambigous 
#    --vep_cache_version s: Cache version, e.g. '90', used to identify cache file.  Optional if not ambiguous
#    --vep_cache_dir s: location of VEP cache directory
#        * if vep_cache_dir is not defined, error.
#        * If vep_cache_dir is a directory, it indicates location of VEP cache 
#        * If vep_cache_dir is a file ending in .tar.gz, will extract its contents into "./vep-cache" and use VEP cache

INPUT_VCF=""  # get output of step 8

if [ -z $VEP_CACHE_DIR ]; then
    echo VEP_CACHE_DIR is not defined.  This step does not work with StrelkaDemo.
    exit 1
fi

ARGS="\
--input_vcf $INPUT_VCF \
--reference_fasta $REFERENCE_FASTA \
--results_dir $RESULTS_DIR \
--assembly $ASSEMBLY \
--vep_cache_version $VEP_CACHE_VERSION \
--vep_cache_dir $VEP_CACHE_DIR \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# Output: results/maf/output.maf
