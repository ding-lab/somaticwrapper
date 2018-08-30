# Testing depth filter using pyvcf's exensible vcf_filter.py framework

DATAD="/data"  # this works if running inside of docker

STRELKA_VCF="$DATAD/origdata/strelka.somatic.snv.all.dbsnp_pass.vcf"
VARSCAN_VCF="$DATAD/dat/varscan.short.vcf"
VARSCAN_INDEL_VCF="$DATAD/origdata/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf"
PINDEL_VCF="$DATAD/dat/pindel.short.vcf"

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
VAF_FILTER_LOCAL="vaf_filter.py"  # filter module
DEPTH_FILTER_LOCAL="depth_filter.py"  # filter module
LENGTH_FILTER_LOCAL="length_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path
CONFIG="--config ../params/vcf_filter_config.ini"

# arguments to depth filter
DEPTH_ARGS="read_depth --debug"

## testing - varscan SNP - OK
CALLER="--caller varscan"
VCF=$VARSCAN_VCF
$MAIN_FILTER --local-script $DEPTH_FILTER_LOCAL $VCF $DEPTH_ARGS $CALLER $CONFIG

### testing - varscan indel - OK
#CALLER="--caller varscan"
#VCF=$VARSCAN_INDEL_VCF
#$MAIN_FILTER --local-script $DEPTH_FILTER_LOCAL $VCF $DEPTH_ARGS $CALLER
#
### testing - strelka - OK
#CALLER="--caller strelka"
#VCF=$STRELKA_VCF
#$MAIN_FILTER --local-script $DEPTH_FILTER_LOCAL $VCF $DEPTH_ARGS $CALLER
#
#### testing - pindel - OK
#CALLER="--caller pindel"
#NAMES="--normal_name pindel.N --tumor_name pindel.T"
#VCF=$PINDEL_VCF
#$MAIN_FILTER --local-script $DEPTH_FILTER_LOCAL $VCF $DEPTH_ARGS $CALLER $NAMES
