# Testing merge filter wrapper
# run_merged_filter.sh exists as a wrapper script so it can easily be called from workflow
source common_config.sh

VCF="StrelkaDemo-results/merged/merged.vcf"

RUN="../../src/vcf_filters/run_merged_filter.sh"
OUT="-"
bash $RUN $VCF $OUT --debug

OUT="tmp/merged1.vcf"
bash $RUN $VCF $OUT --bypass

OUT="tmp/merged2.vcf"
bash $RUN $VCF $OUT --bypass_merge
