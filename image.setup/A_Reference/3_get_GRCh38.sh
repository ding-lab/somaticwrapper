# Download the GRCh38 human reference

OUTD="/image/A_Reference"
mkdir -p $OUTD

echo Saving reference to $OUTD

CWD=`pwd`
cd $OUTD
SRC="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

echo Downloading $SRC
wget $SRC
gunzip -v -f hg38.fa.gz

cd $CWD

echo Preparing reference GRCh38
bash ./prepare_reference.sh $OUTD/hg38.fa
