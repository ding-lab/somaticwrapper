# Download the hg19/GRCh37 human reference
# Not clear this is necessary

OUTD="/image/A_Reference"
mkdir -p $OUTD

echo Saving reference to $OUTD

CWD=`pwd`
cd $OUTD
SRC="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"

echo Downloading $SRC
wget $SRC

echo Downloading twoBitToFa
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod 755 twoBitToFa

./twoBitToFa hg19.2bit hg19.fa
rm hg19.2bit

cd $CWD

echo Preparing reference hg19
bash ./prepare_reference.sh $OUTD/hg19.fa
