# This is run from within a container typically started with ./start_docker.sh

REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=vep_annotate

# In general we want to use the data in the /data directory
#MERGED_VCF="/data/s8_merge_vcf/results/merged/merged.filtered.vcf"
# However, currently we moved merged.filtered.vcf here explicitly for testing
MERGED_VCF="C3N-01649.results/merged/merged.filtered.vcf"

AF_FILTER_CONFIG="/usr/local/somaticwrapper/params/af_filter_config.ini"
CLASS_FILTER_CONFIG="/usr/local/somaticwrapper/params/classification_filter_config.ini"


CACHE_DIR="/image/D_VEP"
CACHE_GZ="/image/D_VEP/vep-cache.90_GRCh37.tar.gz"


# 9 vep_annotate:  
#     --input_vcf s: VCF file to be annotated with vep_annotate.  Required
#     --reference_fasta s: path to reference.  Required
#     --assembly s: either "GRCh37" or "GRCh38", used to identify cache file. Optional if not ambigous 
#     --vep_cache_version s: Cache version, e.g. '90', used to identify cache file.  Optional if not ambiguous
#     --vep_cache_gz: is a file ending in .tar.gz containing VEP cache tarball
#     --vep_cache_dir s: location of VEP cache directory
#         VEP Cache logic:
#         * If vep_cache_dir is defined, it indicates location of VEP cache 
#         * if vep_cache_dir is not defined, and vep_cache_gz is defined, extract vep_cache_gz contents into "./vep-cache" and use VEP cache
#         * if neither vep_cache_dir nor vep_cache_gz defined, will perform online VEP DB lookups
#         NOTE: Online VEP database lookups a) uses online database (so cache isn't installed) b) does not use tmp files
#           It is meant to be used for testing and lightweight applications.  Use the cache for better performance.
#           See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
#     --af_filter_config s: configuration file for af (allele frequency) filter
#     --classification_filter_config s: configuration file for classification filter
#     --bypass: Bypass filter by retaining all reads

ARGS="\
--input_vcf $MERGED_VCF  \
--reference_fasta $REFERENCE_FASTA \
--results_dir $OUTDIR \
--vep_cache_version 90 \
--assembly GRCh37 \
--vep_cache_dir $CACHE_DIR \
--af_filter_config $AF_FILTER_CONFIG \
--classification_filter_config $CLASS_FILTER_CONFIG \
"
#--bypass \

#--vep_cache_gz $CACHE_GZ \
#--vep_output vcf \
#--vep_cache_version 90 \
#--assembly GRCh37 \

# optional:
# --assembly
# --vep_cache_version
# --vep_cache_dir
# --vep_output


BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

