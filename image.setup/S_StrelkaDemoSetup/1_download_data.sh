# Using Strelka test data from https://github.com/Illumina/strelka/tree/master/src/demo/data
# This is a compact dataset which has variants on chromosome "demo20"
# Using NA12892_demo20.bam as Normal
#       NA12891_demo20.bam as Tumor

#github repo: https://github.com/Illumina/strelka
#path of interest: src/demo/data
#
#Downloading of individual folder described here:
#https://stackoverflow.com/questions/7106012/download-a-single-folder-or-directory-from-a-github-repo

# Required
# apt-get install subversion

DATAD_C="/import/StrelkaTestData" # Where test BAMs live
mkdir -p $DATAD_C

SRC="https://github.com/Illumina/strelka"
DIR="src/demo/data"

echo Downloading $DIR from $SRC

svn checkout $SRC/trunk/$DIR/ $DATAD_C

echo Written Strelka test data to $DATAD_C
