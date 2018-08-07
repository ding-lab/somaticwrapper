# Use pyvcf's exensible vcf_flter.py framework



# vcf_filter.py --no-filtered varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf dps --depth-per-sample 100

#VCF="origdata/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
STRELKA_VCF="origdata/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="dat/varscan.short.vcf"
PINDEL_VCF="origdata/pindel.out.current_final.dbsnp_pass.vcf"
PINDEL_VCF="dat/pindel.short.vcf"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
VAF_FILTER_LOCAL="vaf_filter.py"  # filter module
DEPTH_FILTER_LOCAL="depth_filter.py"  # filter module
LENGTH_FILTER_LOCAL="length_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Arguments to VAF filter
SNV_VAF_ARGS="vaf --min_vaf_somatic 0.1 --debug --caller varscan" # --debug"

# arguments to length filter
INDEL_LENGTH_ARGS="--max_length 100 " #--length_debug" 

# pindel depth
DEPTH_ARGS="read_depth --min_depth 200 --caller pindel --debug"
$MAIN_FILTER --local-script $DEPTH_FILTER_LOCAL $PINDEL_VCF $DEPTH_ARGS


# Coverage depth: 20
#$MAIN_FILTER --no-filtered $PINDEL_VCF $DEPTH_ARGS

# Note that depth filter fails for pindel (does not have DP tag)
# need to add new depth filter which works for all callers

# SNV VAF
#$MAIN_FILTER --no-filtered --local-script $VAF_FILTER $VARSCAN_VCF $SNV_VAF_ARGS

# length
#$MAIN_FILTER --no-filtered --local-script $VAF_FILTER $PINDEL_VCF indel_length $INDEL_LENGTH_ARGS $DEPTH_ARGS


# TODO: work through merge logic


