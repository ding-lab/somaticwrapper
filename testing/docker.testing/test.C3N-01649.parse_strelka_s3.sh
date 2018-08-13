# This is run from within a container typically started with ./start_docker.sh

REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"

OUTDIR="./StrelkaDemo.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=3

#    die("Strelka SNV Raw input file not specified \n") unless $strelka_snv_raw;
#    die("strelka_config undefined \n") unless $strelka_config;
#    die("strelka_vcf_filter_config undefined \n") unless $strelka_vcf_filter_config;
#    parse_strelka($results_dir, $job_files_dir, $perl, $gvip_dir, $dbsnp_db, $snpsift_jar, $strelka_snv_raw, $strelka_vcf_filter_config);

# From past TinDaisy run
STRELKA_SNV_RAW="/data/s3_parse_strelka/results/strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_INDEL_VCF="/data/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
VARSCAN_SNV_VCF="/data/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
# from output of parse_pindel
PINDEL_VCF="./StrelkaDemo.results/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf"


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

