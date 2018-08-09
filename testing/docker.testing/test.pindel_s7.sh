# This is run from within a container typically started with ./start_docker.sh

DATAD="../StrelkaDemo.dat"
TUMOR_BAM=$DATAD/StrelkaDemoCase.T.bam
NORMAL_BAM=$DATAD/StrelkaDemoCase.N.bam
REFERENCE_FASTA=$DATAD/demo20.fa
STRELKA_CONFIG=$DATAD/strelka.WES.ini
VARSCAN_CONFIG=$DATAD/varscan.WES.ini
DBSNP_DB=$DATAD/dbsnp-StrelkaDemo.noCOSMIC.vcf.gz
CENTROMERE_BED=$DATAD/ucsc-centromere.GRCh37.bed
#/home/mwyczalk_test/data/docker/data/B_Filter/dbsnp.noCOSMIC.GRCh37.vcf.gz
#/home/mwyczalk_test/data/docker/data
#VEP_CACHE_DIR=/home/mwyczalk_test/data/docker/data/D_VEP
ASSEMBLY="GRCh37"

OUTDIR="./StrelkaDemo.results"
mkdir -p $OUTDIR

SAMPLE="fastq2bam_test"

# Real path is /diskmnt/Projects/Users/hsun/beta_tinDaisy/compare/mgi_sw_C3N-01649/pindel/pindel.out.raw
PINDEL_RAW="/data/pindel/pindel.out.raw"

STEP=7

if [ -z $STEP ]; then
echo Must provide step number.  Usage:
echo    bash run.SomaticWrapper.docker.sh STEP
exit 1
fi

ARGS="\
--reference_fasta $REFERENCE_FASTA \
--pindel_config $PINDEL_CONFIG \
--dbsnp_db $DBSNP_DB \
--no_delete_temp \
--results_dir $OUTDIR \
--pindel_raw $PINDEL_RAW \
"

# final output of step 10 is ./results/merged/merged.vcf (or ./results/vep/output.vcf.vep if --output_vep)

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

