# Download the GRCh38 human reference associated with GDC harmonization
# see https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files



OUTD="/data/image.data/A_Reference"
mkdir -p $OUTD

echo Saving reference to $OUTD

SRC="https://api.gdc.cancer.gov/data/62f23fad-0f24-43fb-8844-990d531947cf"
OUT="GRCh38.d1.vd1.fa.tar.gz"

CWD=`pwd`
cd $OUTD

echo Downloading $SRC
wget $SRC -O $OUT
echo Written to $OUT
tar -zxf $OUT

cd $CWD

echo Preparing reference GRCh38.d1.vd1.fa
bash ./prepare_reference.sh $OUTD/GRCh38.d1.vd1.fa
