source common_config.sh

VEP_VCF="../C3N-01649.test/C3N-01649.results/vep/output.unfiltered.vcf"
AF_CONFIG="../../params/af_filter_config.ini"
CLASS_CONFIG="../../params/classification_filter_config.ini"
#ARGS="--debug --bypass"

RUN="../../src/vcf_filters/run_combined_af_classification_filter.sh"

OUT_VCF="-"

#   bash run_combined_af_classification_filter.sh input.vcf af_config.ini classification_config.ini output.vcf [args]
bash $RUN $VEP_VCF $AF_CONFIG $CLASS_CONFIG $OUT_VCF $ARGS
