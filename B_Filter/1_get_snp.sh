# Download b150 version of all variant VCF, and retain only headers and first 5 columns of VCF

OUTD="/data/B_Filter"
mkdir -p $OUTD

SRC=https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz
OUT=All.5col.vcf.gz

echo downloading $SRC, writing to $OUTD/$OUT

exit

WD=`cwd`

cd $OUTD
wget -q -O - $SRC   | gunzip | perl -lane 'if($_=~/^#/){print $_} else {print join("\t",@F[0..4])}' | bgzip > $OUT

echo Creating TBI file
tabix -p vcf $OUT

cd $WD
