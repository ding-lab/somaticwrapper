# This is run from within a container typically started with ./start_docker.sh

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=3

DBSNP_DB="/image/B_Filter/dbsnp.noCOSMIC.GRCh37.vcf.gz"  # note that this differs from 00-All.brief.pass.cosmic.vcf.gz used in SW.  For present purposes that is OK
STRELKA_SNV_RAW="/data/root/s1_run_strelka/results/strelka/strelka_out/results/passed.somatic.snvs.vcf"
STRELKA_VCF_FILTER_CONFIG="/usr/local/somaticwrapper/params/vcf_filter_config.ini"
STRELKA_CONFIG="../StrelkaDemo.dat/strelka.WES.ini"

ARGS="\
--dbsnp_db $DBSNP_DB \
--strelka_snv_raw $STRELKA_SNV_RAW  \
--results_dir $OUTDIR \
--strelka_config $STRELKA_CONFIG \
--strelka_vcf_filter_config $STRELKA_VCF_FILTER_CONFIG
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

