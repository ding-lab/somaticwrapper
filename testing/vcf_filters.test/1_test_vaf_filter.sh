# Testing vaf filter using pyvcf's exensible vcf_flter.py framework
source common_config.sh

function run_vaf_filter {
CALLER=$1; shift
VCF=$1; shift
XARGS="$@"

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path
FILTER_LOCAL="vaf_filter.py"  # filter module

$MAIN_FILTER --local-script $FILTER_LOCAL $VCF vaf --caller $CALLER $CONFIG $XARGS

}

CONFIG="--config ../../params/vcf_filter_config.ini"

STRELKA_VCF="StrelkaDemo-results/strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf"
#VARSCAN_VCF="$DATAD/root/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
#VARSCAN_INDEL_VCF="$DATAD/root/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
#PINDEL_VCF="$DATAD/root/s7_parse_pindel/results/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf"


run_vaf_filter strelka $STRELKA_VCF --debug

#run_vaf_filter varscan $VARSCAN_VCF --bypass --debug

#run_vaf_filter varscan $VARSCAN_INDEL_VCF --bypass --debug

#run_vaf_filter pindel $PINDEL_VCF --bypass --debug
