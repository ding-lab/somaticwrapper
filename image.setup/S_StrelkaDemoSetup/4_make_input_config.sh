# Create sample input configuration file for Strelka demo analysis

# TODO: this needs to be updated to latest version of config file

DATA_BASE="/data/"   
SAMPLE_NAME="SWtest"
CONFIG_DIR="/data/StrelkaTest/$SAMPLE_NAME/config"

OUT="$CONFIG_DIR/$SAMPLE_NAME.config"  

cat << EOF > $OUT
sample_name = SWtest
reference_fasta = /data/image.data/A_Reference/demo20.fa  
reference_dict = /data/image.data/A_Reference/demo20.dict
tumor_bam = /data/StrelkaTest/data/SWtest/SWtest.T.bam
normal_bam = /data/StrelkaTest/data/SWtest/SWtest.N.bam
assembly = GRCh37
dbsnp_db = /data/image.data/B_Filter/dbsnp-demo.noCOSMIC.vcf.gz

## We will not use VEP cache for testing because it takes a while to install
use_vep_db = 1

# Write annotated data in VEP format (rather than VCF) with gene names 
output_vep = 1

strelka_config = /usr/local/somaticwrapper/config/strelka.WES.ini
varscan_config = /usr/local/somaticwrapper/config/varscan.WES.ini

## command to initiate execution of generated script.  
## 'submit_cmd = cat' will simply print contents of run script, useful for debugging.
## Default value is 'bash'
# submit_cmd = cat

# VEP-annotate intermediate output files (testing)
annotate_intermediate = 1

## MGI development requires a non-default location for somaticwrapper
# sw_dir = /gscuser/mwyczalk/projects/SomaticWrapper/somaticwrapper
EOF

echo Written configuration to $OUT
