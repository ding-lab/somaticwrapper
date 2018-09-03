# This is run from within a container typically started with ./start_docker.sh
# Testing processing of chrom 22, which we know has a variant after filtering

RESULTS_DIR="./C3N-01649.results"
mkdir -p $RESULTS_DIR

SAMPLE="C3N-01649"

STEP="run_pindel"

NORMAL_BAM="/GDC_import/data/9c91aa54-a11b-4316-a8a9-3aa38aa09771/CPT0088680009.WholeExome.RP-1303.bam"
TUMOR_BAM="/GDC_import/data/c881268c-c1e5-4803-ba10-64f03f3fdced/CPT0088640009.WholeExome.RP-1303.bam"
REFERENCE_FASTA="/image/A_Reference/Homo_sapiens_assembly19.fasta"

CENTROMERE_BED="/usr/local/somaticwrapper/testing/StrelkaDemo.dat/ucsc-centromere.GRCh37.bed"
PINDEL_CONFIG="/usr/local/src/params/pindel.WES.ini"

#5 run_pindel:
#    --tumor_bam s:  path to tumor BAM.  Required
#    --normal_bam s: path to normal BAM.  Required
#    --reference_fasta s: path to reference.  Required
#    --pindel_config s: path to pindel.ini file.  Required
#    --no_delete_temp : if defined, do not delete temp files
#    --centromere_bed s: path to BED file describing centromere regions to exclude for pindel analysis. Optional
#    --pindel_chrom c: defines which chrom to process; default is to process all of them.  Optional

ARGS="\
--tumor_bam $TUMOR_BAM \
--normal_bam $NORMAL_BAM \
--reference_fasta $REFERENCE_FASTA \
--pindel_config $PINDEL_CONFIG \
--centromere_bed $CENTROMERE_BED \
--results_dir $RESULTS_DIR \
--no_delete_temp \
--pindel_chrom 22 \
"  

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

# output: results/pindel/pindel_out/pindel-raw.dat
