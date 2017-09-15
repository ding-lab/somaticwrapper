
DATD="/data/A_Reference"
DAT="GRCh37-lite.fa"

CWD=`pwd`
cd $DATD

samtools faidx $DAT
/usr/bin/bwa index $DAT

cd $CWD
