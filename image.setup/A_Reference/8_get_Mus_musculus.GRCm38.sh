# Download the Mouse reference
SRC="ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
FAZ="Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
FA="Mus_musculus.GRCm38.dna_sm.primary_assembly.fa"

OUTD="/data/A_Reference"

echo Saving reference to $OUTD

CWD=`pwd`
cd $OUTD

wget $SRC
gunzip -v -f $FAZ

cd $CWD

echo Preparing reference $FA
bash ./prepare_reference.sh $OUTD/$FA
