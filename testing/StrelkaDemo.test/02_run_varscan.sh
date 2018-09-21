DATAD="/usr/local/somaticwrapper/testing/StrelkaDemo.dat"
source project_config.sh $DATAD

STEP="run_varscan"

#2 run_varscan:
#    --tumor_bam s:  path to tumor BAM.  Required
#    --normal_bam s: path to normal BAM.  Required
#    --reference_fasta s: path to reference.  Required
#    --varscan_config s: path to varscan.ini file.  Required
#    --results_dir s: Per-sample analysis results location. Often same as sample name [.] 

ARGS=" \
--tumor_bam $TUMOR_BAM \
--normal_bam $NORMAL_BAM \
--reference_fasta $REFERENCE_FASTA \
--varscan_config $VARSCAN_CONFIG \
--results_dir $RESULTS_DIR \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

#Final results: SNV: results/varscan/varscan_out/varscan.out.som_snv.vcf
#             INDEL: results/varscan/varscan_out/varscan.out.som_indel.vcf

