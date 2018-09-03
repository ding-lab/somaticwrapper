# This is run from within a container typically started with ./start_docker.sh

REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"
PINDEL_CONFIG="../StrelkaDemo.dat/pindel.WES.ini"
DBSNP_DB="/image/B_Filter/dbsnp.noCOSMIC.GRCh37.vcf.gz"  # note that this differs from 00-All.brief.pass.cosmic.vcf.gz used in SW.  For present purposes that is OK
PINDEL_RAW="/data/root/s5_run_pindel/results/pindel/pindel_out/pindel-raw.dat"

# for testing, renamed pindel.N, pindel.T -> NORMAL, TUMOR
# future pindel runs will use this naming convention
PINDEL_RAW="/usr/local/somaticwrapper/testing/C3N-01649.test/C3N-01649.results/pindel/testdat/pindel-raw.renamed.dat"
PINDEL_VCF_FILTER_CONFIG="/usr/local/somaticwrapper/params/pindel-vcf_filter_config.ini"

OUTDIR="./C3N-01649.results"
mkdir -p $OUTDIR

SAMPLE="C3N-01649"

STEP=7

ARGS="\
--reference_fasta $REFERENCE_FASTA \
--pindel_config $PINDEL_CONFIG \
--dbsnp_db $DBSNP_DB \
--results_dir $OUTDIR \
--pindel_raw $PINDEL_RAW \
--pindel_vcf_filter_config $PINDEL_VCF_FILTER_CONFIG \
--no_delete_temp \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

