#!/bin/bash

# Usage: download_cosmic.sh REF VER
# e.g., 
# REF="grch37"
# VER="v82"
REF=$1
VER=$2


# Download data from COSMIC
# You need to be registered to access the data: http://cancer.sanger.ac.uk/cosmic
# Username and password are stored in COSMIC_credentials.dat
#   If the file does not exist, just run this script to create it
# Saves to ./dat

DAT="./COSMIC_credentials.dat"

OUTD="/data/image.data/B_Filter"
mkdir -p $OUTD

if [ ! -f $DAT ]; then
cat <<EOF > $DAT
# Use username/password you registered at http://cancer.sanger.ac.uk/cosmic
export COSMIC_USERNAME="cosmic_username"
export SSHPASS="cosmic_password"
EOF

echo Created template file $DAT 
echo Please edit this file to add your credentials and run this script again
exit

else
source ./COSMIC_credentials.dat
fi

# needs sshpass for automated download (apt-get install sshpass)
#
# An alternative is to simply download it by hand like,
#   sftp "$COSMIC_USERNAME"@sftp-cancer.sanger.ac.uk
#   <enter your password>
#   cd /cosmic/$REF/cosmic/$VER/VCF
#   get CosmicCodingMuts.vcf.gz

COSD="/cosmic/$REF/cosmic/$VER/VCF"
VCFGZ="CosmicCodingMuts.vcf.gz"             # file as downloaded
VCFBGZ="CosmicCodingMuts.$REF.$VER.vcf.gz"  # file as processed
VCF="$COSD/$VCFGZ"

if [ -f $VCFGZ ]; then
echo File $VCFGZ exists.  Please delete if you want to re-download
exit
fi

if [ -f $VCFBGZ ]; then
echo File $VCFBGZ exists.  Please delete if you want to re-download
exit
fi

echo Downloading $VCF to $OUTD:
echo sftp "$COSMIC_USERNAME"@sftp-cancer.sanger.ac.uk
# if the download below doesn't work, run sftp by hand.  May need to verify authenticity of host first time doing this

sshpass -e sftp -oBatchMode=no -b - "$COSMIC_USERNAME"@sftp-cancer.sanger.ac.uk << EOF
lcd $OUTD
get $VCF
bye
EOF

# Testing to see if downloaded file exists.  If not, provide advice
if [ ! -f $OUTD/$VCFGZ ]; then

cat <<EOF
Error downloading data.  Please confirm login credentials by running the command,
    sftp "$COSMIC_USERNAME"@sftp-cancer.sanger.ac.uk

In some cases need to verify authenticity of host just the first time this command is run.
EOF
exit

else

echo $OUTD/$VCFGZ downloaded successfully

fi

# Now convert to bgz format and index.  Note that extension needs to be .gz for downstream
# java snpsift code to work
gunzip -c $OUTD/$VCFGZ | bgzip > $OUTD/$VCFBGZ
tabix -p vcf $OUTD/$VCFBGZ

echo Written to $OUTD/$VCFBGZ

rm $OUTD/$VCFGZ

