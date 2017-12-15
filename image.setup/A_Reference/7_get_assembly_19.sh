# This reference looks like this in the BAM header:
# UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta	AS:GRCh37
# Chrom 1 is labeled as '1'

OUTD="/data/image.data/A_Reference"
mkdir -p $OUTD

echo Saving reference to $OUTD

function get_bundle {
F=$1

SRC="http://archive.broadinstitute.org/ftp/pub/seq/references/$F"

echo Downloading $SRC
wget $SRC
#gunzip -v -f $F

}

CWD=`pwd`
cd $OUTD

get_bundle Homo_sapiens_assembly19.fasta
get_bundle Homo_sapiens_assembly19.fasta.fai
get_bundle Homo_sapiens_assembly19.dict

cd $CWD

