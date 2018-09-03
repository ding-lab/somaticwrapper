# This is run from within a container typically started with ./start_docker.sh

REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=vcf_2_maf

# In general we want to use the data in the /data directory
#MERGED_VCF="/data/s8_merge_vcf/results/merged/merged.filtered.vcf"
# However, currently we moved merged.filtered.vcf here explicitly for testing
MERGED_VCF="./C3N-01649.results/vep/output.vcf"

AF_FILTER_CONFIG="/usr/local/somaticwrapper/params/af_filter_config.ini"
CLASS_FILTER_CONFIG="/usr/local/somaticwrapper/params/classification_filter_config.ini"


CACHE_DIR="/image/D_VEP"
CACHE_GZ="/image/D_VEP/vep-cache.90_GRCh37.tar.gz"

#10 vcf_2_maf:
#    --input_vcf s: VCF file to be annotated with vep_annotate.  Required
#    --reference_fasta s: path to reference.  Required
#    --assembly s: either "GRCh37" or "GRCh38", used to identify cache file. Optional if not ambigous 
#    --vep_cache_version s: Cache version, e.g. '90', used to identify cache file.  Optional if not ambiguous
#    --vep_cache_gz: is a file ending in .tar.gz containing VEP cache tarball
#    --vep_cache_dir s: location of VEP cache directory
#        VEP Cache logic:
#        * If vep_cache_dir is defined, it indicates location of VEP cache 
#        * if vep_cache_dir is not defined, and vep_cache_gz is defined, extract vep_cache_gz contents into "./vep-cache" and use VEP cache
#        * if neither vep_cache_dir nor vep_cache_gz defined, error.  vcf_2_maf does not support online vep_cache lookups
#    --exac:  ExAC database to pass as --f_exac for annotation
#
#    --results_dir s: Per-sample analysis results location. Often same as sample name [.] 
#    --skip s: If defined, skip this step and print argument as reason for skipping.  Helpful for interaction with CWL workflow.

ARGS="\
--input_vcf $MERGED_VCF  \
--reference_fasta $REFERENCE_FASTA \
--results_dir $OUTDIR \
--vep_cache_version 90 \
--assembly GRCh37 \
--vep_cache_dir $CACHE_DIR \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP
