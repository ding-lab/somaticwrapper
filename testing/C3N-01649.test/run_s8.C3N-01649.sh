# This is run from within a container typically started with ./start_docker.sh

REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=8

DATAD="/data/root"
STRELKA_SNV_VCF="$DATAD/s3_parse_strelka/results/strelka/filter_out/dbsnp_pass.filtered.vcf"
VARSCAN_SNV_VCF="$DATAD/s4_parse_varscan/results/varscan/filter_out/dbsnp_pass.filtered.vcf"
VARSCAN_INDEL_VCF="$DATAD/s4_parse_varscan/results/varscan/filter_out/dbsnp_pass.filtered.vcf"
PINDEL_VCF="$DATAD/s7_parse_pindel/results/pindel/filter_out/dbsnp_pass.filtered.vcf"

ARGS="\
--reference_fasta $REFERENCE_FASTA \
--strelka_snv_vcf $STRELKA_SNV_VCF  \
--varscan_indel_vcf $VARSCAN_INDEL_VCF  \
--varscan_snv_vcf $VARSCAN_SNV_VCF  \
--pindel_vcf $PINDEL_VCF \
--results_dir $OUTDIR \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

#   1554 varscan
#    245 pindel
#     88 varindel
#     24 strelka
#      9 strelka-varscan

