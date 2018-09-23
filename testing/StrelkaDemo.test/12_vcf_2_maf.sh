DATAD="/usr/local/somaticwrapper/testing/StrelkaDemo.dat"
source project_config.sh $DATAD

# This step is not implemented for StrelkaDemo because 
# vcf2maf.pl requires a VEP cache, which is in general not installed for the demo
# With a cache available (either installed or as .tar.gz) 

STEP="vcf_2_maf"

# vcf_2_maf:
#     --input_vcf s: VCF file to be annotated with vep_annotate.  Required
#     --reference_fasta s: path to reference.  Required
#     --assembly s: either "GRCh37" or "GRCh38", used to identify cache file. Optional if not ambigous 
#     --vep_cache_version s: Cache version, e.g. '90', used to identify cache file.  Optional if not ambiguous
#     --vep_cache_gz: is a file ending in .tar.gz containing VEP cache tarball
#     --vep_cache_dir s: location of VEP cache directory
#         VEP Cache logic:
#         * If vep_cache_dir is defined, it indicates location of VEP cache 
#         * if vep_cache_dir is not defined, and vep_cache_gz is defined, extract vep_cache_gz contents into "./vep-cache" and use VEP cache
#         * if neither vep_cache_dir nor vep_cache_gz defined, error.  vcf_2_maf does not support online vep_cache lookups
#     --exac:  ExAC database to pass to vcf_2_maf.pl as --filter-vcf for custom annotation

INPUT_VCF="results/merged/merged.filtered.vcf"  

if [ -z $VEP_CACHE_DIR ]; then
    echo VEP_CACHE_DIR is not defined.  This step does not work with StrelkaDemo.
    # Getopt::Long does not seem to handle spaces in arguments well.  Passing awkward error message
    SKIP="--skip VEP_CACHE.not.defined"
fi

ARGS="\
--input_vcf $INPUT_VCF \
--reference_fasta $REFERENCE_FASTA \
--results_dir $RESULTS_DIR \
--assembly $ASSEMBLY \
--vep_cache_version $VEP_CACHE_VERSION \
--vep_cache_dir $VEP_CACHE_DIR \
$SKIP \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# Output: results/maf/output.maf
