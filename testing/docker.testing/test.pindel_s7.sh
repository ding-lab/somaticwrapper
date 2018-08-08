# This is run from within a container typically started with ./start_docker.sh

DATAD="../StrelkaDemo.dat"
TUMOR_BAM=$DATAD/StrelkaDemoCase.T.bam
NORMAL_BAM=$DATAD/StrelkaDemoCase.N.bam
REFERENCE_FASTA=$DATAD/demo20.fa
STRELKA_CONFIG=$DATAD/strelka.WES.ini
VARSCAN_CONFIG=$DATAD/varscan.WES.ini
PINDEL_CONFIG=$DATAD/pindel.WES.ini
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
--tumor_bam $TUMOR_BAM \
--normal_bam $NORMAL_BAM \
--reference_fasta $REFERENCE_FASTA \
--strelka_config $STRELKA_CONFIG \
--varscan_config $VARSCAN_CONFIG \
--pindel_config $PINDEL_CONFIG \
--dbsnp_db $DBSNP_DB \
--nooutput_vep \
--assembly $ASSEMBLY \
--is_strelka2  \
--no_delete_temp \
--centromere_bed $CENTROMERE_BED \
--results_dir $OUTDIR \
--strelka_snv_raw $OUTDIR/strelka/strelka_out/results/variants/somatic.snvs.vcf.gz \
--strelka_snv_vcf $OUTDIR/strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf \
--varscan_indel_raw $OUTDIR/varscan/varscan_out/varscan.out.som_indel.vcf \
--varscan_snv_raw $OUTDIR/varscan/varscan_out/varscan.out.som_snv.vcf \
--varscan_indel_vcf $OUTDIR/varscan/filter_out/varscan.out.som_indel.Somatic.hc.vcf \
--varscan_snv_vcf $OUTDIR/varscan/filter_out/varscan.out.som_snv.Somatic.hc.vcf \
--pindel_raw $PINDEL_RAW \
--pindel_vcf $OUTDIR/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf \
--input_vcf $OUTDIR/merged/merged.vcf \
"

# final output of step 10 is ./results/merged/merged.vcf (or ./results/vep/output.vcf.vep if --output_vep)

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

