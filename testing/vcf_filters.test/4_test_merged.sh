# Testing merge filter wrapper
# run_merged_filter.sh exists as a wrapper script so it can easily be called from workflow
source common_config.sh

VCF="$DATAD/root/s8_merge_vcf/results/merged/merged.vcf"

OUT="-"


RUN="../../src/vcf_filters/run_merged_filter.sh"
bash $RUN $VCF $OUT --debug

