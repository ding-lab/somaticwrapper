DATAD="/usr/local/somaticwrapper/testing/StrelkaDemo.dat"
source project_config.sh $DATAD

STEP="run_pindel"

#5 run_pindel:
#    --tumor_bam s:  path to tumor BAM.  Required
#    --normal_bam s: path to normal BAM.  Required
#    --reference_fasta s: path to reference.  Required
#    --pindel_config s: path to pindel.ini file.  Required
#    --results_dir s: Per-sample analysis results location. Often same as sample name [.] 
#    --no_delete_temp : if defined, do not delete temp files
#    --centromere_bed s: path to BED file describing centromere regions to exclude for pindel analysis.  

ARGS="\
--tumor_bam $TUMOR_BAM \
--normal_bam $NORMAL_BAM \
--reference_fasta $REFERENCE_FASTA \
--centromere_bed $CENTROMERE_BED \
--results_dir $RESULTS_DIR \
"  
#--no_delete_temp \

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# output: results/pindel/pindel_out/pindel-raw.dat
