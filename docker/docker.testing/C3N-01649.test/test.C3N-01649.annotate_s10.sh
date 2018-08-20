# This is run from within a container typically started with ./start_docker.sh

REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=10

MERGED_VCF="$OUTDIR/merged/merged.filtered.vcf"

CACHE_DIR="/image/D_VEP"
#CACHE_DIR="/image/D_VEP/vep-cache.90_GRCh37.tar.gz"

ARGS="\
--reference_fasta $REFERENCE_FASTA \
--input_vcf $MERGED_VCF  \
--results_dir $OUTDIR \
--vep_cache_dir $CACHE_DIR \
--vep_output vcf \
--vep_cache_version 90 \
--assembly GRCh37 \
"
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

