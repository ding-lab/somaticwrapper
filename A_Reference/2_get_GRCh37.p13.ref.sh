# Download the GRCh37 human reference
# According to Cyriac (https://www.biostars.org/p/119295/#119308), old GRCh37 uses the '1' format, but newer versions (as well as hg19) use the 'chr1'.
# So, we will download the newer GRCh37 from here: ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/seqs_for_alignment_pipelines/


OUTD="/data/A_Reference"
mkdir -p $OUTD

echo Saving reference to $OUTD

CWD=`pwd`
cd $OUTD
SRC="ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37.p13/seqs_for_alignment_pipelines/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz"

echo Downloading $SRC
wget $SRC
gunzip -v -f GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna.gz
mv GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna GRCh37.p13.fa

cd $CWD

echo Preparing reference GRCh37.p13.fa
bash ./prepare_reference.sh $OUTD/GRCh37.p13.fa
