# Create sample input configuration file for Strelka demo analysis

DATA_BASE="/data"   
SAMPLE_NAME="SWtest"
SAMPLE_DIR="$DATA_BASE/data/$SAMPLE_NAME"

OUT="$SAMPLE_DIR/sw.config"  # /data/data/SWtest/sw.config

cat << EOF > $OUT
sample_name = SWtest
reference_fasta = /data/A_Reference/demo20.fa  
tumor_bam = /data/data/SWtest/SWtest.T.bam
normal_bam = /data/data/SWtest/SWtest.N.bam

# MGI development requires a non-default location for somaticwrapper
sw_dir = /gscuser/mwyczalk/projects/SomaticWrapper/somaticwrapper
EOF

echo Written configuration to $OUT
