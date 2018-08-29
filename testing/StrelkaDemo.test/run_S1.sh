source project_config.sh

STEP="run_strelka"

#1 run_strelka:
#    --tumor_bam s:  path to tumor BAM.  Required
#    --normal_bam s: path to normal BAM.  Required
#    --reference_fasta s: path to reference.  Required
#    --strelka_config s: path to strelka.ini file.  Required
#    --results_dir s: Per-sample analysis results location. Often same as sample name [.] 
#    --is_strelka2: run stelka2 instead of strelka version 1 

ARGS="\
--tumor_bam $TUMOR_BAM \
--normal_bam $NORMAL_BAM \
--reference_fasta $REFERENCE_FASTA \
--strelka_config $STRELKA_CONFIG \
--results_dir $RESULTS_DIR \
--is_strelka2 \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# Output (--is_strelka2): results/strelka/strelka_out/results/variants/somatic.snvs.vcf.gz
# With Strelka1: results/strelka/strelka_out/results/passed.somatic.snvs.vcf
