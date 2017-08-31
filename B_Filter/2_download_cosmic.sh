# Download data from COSMIC
# You need to be registered to access the data: http://cancer.sanger.ac.uk/cosmic
# Username and password are stored in COSMIC_credentials.dat
#   If the file does not exist, just run this script to create it
# Saves to ./dat

DAT="./COSMIC_credentials.dat"

OUTD="/data/B_Filter"
mkdir -p $OUTD

if [ ! -f $DAT ]; then
cat <<EOF > $DAT
# Use username/password you registered at http://cancer.sanger.ac.uk/cosmic
COSMIC_USERNAME="cosmic_username"
SSHPASS="cosmic_password"
EOF

echo Please edit the file $DAT with your credentials and run this script again
exit

else
source ./COSMIC_credentials.dat
fi

echo $COSMIC_USERNAME

REF="grch37"

# Cosmic Version
VER="v82"

# needs sshpass for automated download (apt-get install sshpass)
#
# An alternative is to simply download it by hand like,
#   sftp "$COSMIC_USERNAME"@sftp-cancer.sanger.ac.uk
#   <enter your password>
#   cd /cosmic/$REF/cosmic/$VER/VCF
#   get CosmicCodingMuts.vcf.gz

VCF="/cosmic/$REF/cosmic/$VER/VCF/CosmicCodingMuts.vcf.gz"
echo Downloading $VCF to $DAT

export SSHPASS="SauCer+7067"
sshpass -e sftp -oBatchMode=no -b - "$COSMIC_USERNAME"@sftp-cancer.sanger.ac.uk << EOF

lcd $OUTD
get $VCF

bye
EOF

