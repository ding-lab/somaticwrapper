# Usage:
# bash get_snp.sh SRC OUT
# Where for example,
# SRC=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-All.vcf.gz
# OUT=human_9606_b150_GRCh37p13.All.5col.vcf.gz

SRC=$1
OUT=$2



OUTD="/data/image.data/B_Filter"
mkdir -p $OUTD

# See here for file types and sources of links
# https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/#table-1


# GRCh37 common
#SRC=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-common_all.vcf.gz
#OUT=human_9606_b150_GRCh37p13.Common.5col.vcf.gz
# This has '1' notation 

if [ -e $OUTD/$OUT ]; then

echo File $OUTD/$OUT exists.  Please delete it if you want to download again.
exit

fi


echo downloading $SRC, writing to $OUTD/$OUT

WD=`pwd`
cd $OUTD

wget -q -O - $SRC   | gunzip | perl -lane 'if($_=~/^#/){print $_} else {print join("\t",@F[0..4])}' | bgzip > $OUT

echo Creating TBI file
tabix -p vcf $OUT

cd $WD
