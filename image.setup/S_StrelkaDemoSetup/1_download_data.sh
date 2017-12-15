
#github repo: https://github.com/Illumina/strelka
#path of interest: src/demo/data
#
#Downloading of individual folder described here:
#https://stackoverflow.com/questions/7106012/download-a-single-folder-or-directory-from-a-github-repo

# Requied
# apt-get install subversion

OUTD="/data/image.data/S_StrelkaTestData"
mkdir -p $OUTD

SRC="https://github.com/Illumina/strelka"
DIR="src/demo/data"

echo Downloading $DIR from $SRC

svn checkout $SRC/trunk/$DIR/ $OUTD

echo Written Strelka test data to $OUTD
