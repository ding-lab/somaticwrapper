#!/bin/bash

REF="grch38"
VER="v82"

# needs sshpass for automated download (apt-get install sshpass)
#
# An alternative is to simply download it by hand like,
#   sftp "$COSMIC_USERNAME"@sftp-cancer.sanger.ac.uk
#   <enter your password>
#   cd /cosmic/$REF/cosmic/$VER/VCF
#   get CosmicCodingMuts.vcf.gz

bash ./download_cosmic.sh $REF $VER

