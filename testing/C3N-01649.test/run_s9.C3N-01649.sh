# This is run from within a container typically started with ./start_docker.sh

REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=9

#MERGED_VCF="$OUTDIR/merged/merged.filtered.vcf"
MERGED_VCF="/data/s8_merge_vcf/results/merged/merged.filtered.vcf"

CACHE_DIR="/image/D_VEP"
CACHE_GZ="/image/D_VEP/vep-cache.90_GRCh37.tar.gz"

# TODO: test annotation
# * cache-dir
# * cache-gz
# * VEP DB

# Also test passing of af_gnomad and af_exac

#9 vep_annotate:  
#    --input_vcf s: VCF file to be annotated with vep_annotate.  Required
#    --reference_fasta s: path to reference.  Required
#    --results_dir s: Per-sample analysis results location. Often same as sample name [.] 
#    --assembly s: either "GRCh37" or "GRCh38", used to identify cache file. Optional if not ambigous 
#    --vep_cache_version s: Cache version, e.g. '90', used to identify cache file.  Optional if not ambiguous
#    --vep_cache_gz: is a file ending in .tar.gz containing VEP cache tarball
#    --vep_cache_dir s: location of VEP cache directory
#        VEP Cache logic:
#        * If vep_cache_dir is defined, it indicates location of VEP cache 
#        * if vep_cache_dir is not defined, and vep_cache_gz is defined, extract vep_cache_gz contents into "./vep-cache" and use VEP cache
#        * if neither vep_cache_dir nor vep_cache_gz defined, will perform online VEP DB lookups
#        NOTE: Online VEP database lookups a) uses online database (so cache isn't installed) b) does not use tmp files
#          It is meant to be used for testing and lightweight applications.  Use the cache for better performance.
#          See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
#    --vep_output: Define output format after annotation.  Allowed values: vcf, vep.  [vcf]

ARGS="\
--input_vcf $MERGED_VCF  \
--reference_fasta $REFERENCE_FASTA \
--results_dir $OUTDIR \
--vep_cache_version 90 \
--assembly GRCh37 \
--vep_cache_dir $CACHE_DIR \
"
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

#   1554 varscan
#    245 pindel
#     88 varindel
#     24 strelka
#      9 strelka-varscan

