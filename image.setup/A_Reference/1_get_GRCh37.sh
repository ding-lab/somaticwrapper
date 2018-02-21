# Download the GRCh37 human reference 

OUTD="/image/A_Reference"
mkdir -p $OUTD

echo Saving reference to $OUTD

CWD=`pwd`
cd $OUTD

wget ftp://genome.wustl.edu/pub/reference/GRCh37-lite/GRCh37-lite.fa.gz
gunzip -v -f GRCh37-lite.fa.gz

cd $CWD

echo Preparing reference GRCh37-lite.fa
bash ./prepare_reference.sh $OUTD/GRCh37-lite.fa
