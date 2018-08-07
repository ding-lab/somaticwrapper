# Testing pyvcf's exensible vcf_flter.py framework

DATAD="/data"  # this works if running inside of docker

STRELKA_VCF="$DATAD/origdata/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="$DATAD/dat/varscan.short.vcf"
VARSCAN_INDEL_VCF="$DATAD/origdata/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/dat/pindel.short.vcf"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
VAF_FILTER_LOCAL="vaf_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path

# Arguments to VAF filter
SNV_VAF_ARGS="vaf --min_vaf_somatic 0.1 --debug" # --debug"

## testing - varscan SNP - OK
#CALLER="--caller varscan"
#$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $VARSCAN_VCF $SNV_VAF_ARGS $CALLER

# testing - varscan INDEL - OK
#CALLER="--caller varscan"
#$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $VARSCAN_INDEL_VCF $SNV_VAF_ARGS $CALLER

## testing - pindel.  Note that pindel has different sample names for now
#CALLER="--caller pindel" - OK
#NAMES="--normal_name pindel.N --tumor_name pindel.T"
#$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $PINDEL_VCF $SNV_VAF_ARGS $CALLER $NAMES

# testing - strelka - OK
CALLER="--caller strelka"
$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $STRELKA_VCF $SNV_VAF_ARGS $CALLER

