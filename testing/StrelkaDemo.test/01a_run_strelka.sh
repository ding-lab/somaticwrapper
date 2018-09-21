DATAD="/usr/local/somaticwrapper/testing/StrelkaDemo.dat"
source project_config.sh $DATAD

STEP="run_strelka"

# run_strelka:
#     --tumor_bam s:  path to tumor BAM.  Required
#     --normal_bam s: path to normal BAM.  Required
#     --reference_fasta s: path to reference.  Required
#     --strelka_config s: path to strelka.ini file.  Required

ARGS="\
--tumor_bam $TUMOR_BAM \
--normal_bam $NORMAL_BAM \
--reference_fasta $REFERENCE_FASTA \
--strelka_config $STRELKA_CONFIG \
--results_dir $RESULTS_DIR \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# Output (--is_strelka2): results/strelka/strelka_out/results/variants/somatic.snvs.vcf.gz
# With Strelka1: results/strelka/strelka_out/results/passed.somatic.snvs.vcf
