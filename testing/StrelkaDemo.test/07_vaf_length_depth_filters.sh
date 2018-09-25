DATAD="/usr/local/somaticwrapper/testing/StrelkaDemo.dat"
source project_config.sh $DATAD

STEP="vaf_length_depth_filters"

# vaf_length_depth_filters: apply VAF, indel length, and read depth filters to a VCF
#     --input_vcf s: VCF file to process.  Required
#     --output_vcf s: Name of output VCF file (written to results_dir/vaf_length_depth_filters/output_vcf).  Required.
#     --caller s: one of strelka, pindel, varscan - ignoring this, will come from config file
#     --vcf_filter_config: Configuration file for VCF filtering (depth, VAF, read count).  Required
#     --bypass_vaf: skip VAF filter
#     --bypass_length: skip length filter
#     --bypass_depth: skip depth filter
#     --bypass: skip all filters
#     --debug: print out processing details to STDERR

function run_vld_filter {
# Usage: run_vld_filter INPUT_VCF OUTPUT_VCF XARG

ARGS="\
--input_vcf $1 \
--output_vcf $2 \
--vcf_filter_config $3 \
--results_dir $RESULTS_DIR \
$5 \
"

BIN="/usr/local/somaticwrapper/SomaticWrapper.pl"
perl $BIN $ARGS $STEP

}

INPUT_VCF="results/strelka2/strelka_out/results/variants/somatic.snvs.vcf.gz"
run_vld_filter $INPUT_VCF strelka.snv.vcf $STRELKA_VCF_FILTER_CONFIG

INPUT_VCF="results/varscan/varscan_out/varscan.out.som_snv.vcf"
run_vld_filter $INPUT_VCF varscan.snv.vcf $VARSCAN_VCF_FILTER_CONFIG

INPUT_VCF="results/varscan/varscan_out/varscan.out.som_indel.vcf" 
run_vld_filter $INPUT_VCF varscan.indel.vcf $VARSCAN_VCF_FILTER_CONFIG --debug

INPUT_VCF="results/pindel/filter_out/pindel-raw.dat.CvgVafStrand_pass.Homopolymer_pass.vcf"
run_vld_filter $INPUT_VCF pindel.indel.vcf $PINDEL_VCF_FILTER_CONFIG --bypass

