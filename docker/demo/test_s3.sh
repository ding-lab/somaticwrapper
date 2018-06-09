export TUMOR_BAM=/data/StrelkaDemoCase.T.bam
export NORMAL_BAM=/data/StrelkaDemoCase.N.bam
export REFERENCE_FASTA=/data/demo20.fa
export STRELKA_CONFIG=/data/strelka.WES.ini
export PINDEL_CONFIG=/data/pindel.WES.ini
export DBSNP_DB=/data/dbsnp-StrelkaDemo.noCOSMIC.vcf.gz
export CENTROMERE_BED=/data/ucsc-centromere.GRCh37.bed
#export VEP_CACHE_DIR=/home/mwyczalk_test/data/docker/data/D_VEP

export OUTPUT_DIR=./results  #/diskmnt/Projects/TinDaisy/results

export VARSCAN_CONFIG=../../params/varscan.WES.ini

SNV_RAW="/usr/local/somaticwrapper/docker/demo/results/strelka/strelka_out/results/variants/somatic.snvs.vcf.gz"

# try to have all output go to output_dir
mkdir -p $OUTPUT_DIR

ARGS="--strelka_snv_raw $SNV_RAW \
--dbsnp_db $DBSNP_DB \
--strelka_config $STRELKA_CONFIG "

perl ../../SomaticWrapper.pl $ARGS --results_dir $OUTPUT_DIR 3

# run_cwl_S4.sh
#$RABIX $RABIX_ARGS $CWL -- " \
#--varscan_indel_raw /Users/mwyczalk/Projects/Rabix/SomaticWrapper.CWL1/results/s2_run_varscan-2018-03-24-133230.292/root/varscan/varscan_out/varscan.out.som_indel.vcf \
#--varscan_snv_raw /Users/mwyczalk/Projects/Rabix/SomaticWrapper.CWL1/results/s2_run_varscan-2018-03-24-133230.292/root/varscan/varscan_out/varscan.out.som_snv.vcf "
