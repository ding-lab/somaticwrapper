# Download "UCSC_hg19" which comes from Broad Institute bundle.
# It is distinguished by having `chrM` as the first chrom.
# This reference used for 01BR001 test case
# We'll download the fai and dict files so there's no need to prepare the reference ourselves

#OUTD="/data/A_Reference"
OUTD="/Users/mwyczalk/src/SomaticWrapper/data//A_Reference"
mkdir -p $OUTD

echo Saving reference to $OUTD

function get_bundle {
F=$1

SRC="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/$F"

echo Downloading $SRC
wget $SRC
gunzip -v -f $F

}

CWD=`pwd`
cd $OUTD

get_bundle ucsc.hg19.dict.gz
get_bundle ucsc.hg19.fasta.fai.gz
get_bundle ucsc.hg19.fasta.gz

cd $CWD

