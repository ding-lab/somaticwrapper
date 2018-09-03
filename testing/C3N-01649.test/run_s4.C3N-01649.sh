# This is run from within a container typically started with ./start_docker.sh

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=4

DBSNP_DB="/image/B_Filter/dbsnp.noCOSMIC.GRCh37.vcf.gz"  # note that this differs from 00-All.brief.pass.cosmic.vcf.gz used in SW.  For present purposes that is OK
VARSCAN_SNV_RAW="/data/root/s2_run_varscan/results/varscan/varscan_out/varscan.out.som_snv.vcf"
VARSCAN_INDEL_RAW="/data/root/s2_run_varscan/results/varscan/varscan_out/varscan.out.som_indel.vcf"
VARSCAN_VCF_FILTER_CONFIG="/usr/local/somaticwrapper/params/vcf_filter_config.ini"
VARSCAN_CONFIG="../StrelkaDemo.dat/varscan.WES.ini"

ARGS="\
--dbsnp_db $DBSNP_DB \
--varscan_snv_raw $VARSCAN_SNV_RAW  \
--varscan_indel_raw $VARSCAN_INDEL_RAW  \
--results_dir $OUTDIR \
--varscan_config $VARSCAN_CONFIG \
--varscan_vcf_filter_config $VARSCAN_VCF_FILTER_CONFIG \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

