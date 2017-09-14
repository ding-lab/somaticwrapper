
DATD="/data/A_Reference"
DAT="GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna"

CWD=`pwd`
cd $DATD

samtools faidx $DAT
/usr/bin/bwa index $DAT

cd $CWD
