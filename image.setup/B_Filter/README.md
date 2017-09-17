# Create list of common SNPs

1. Download dbsnp database
2. Download COSMIC database
3. Remove from dbsnp database those mutations which exist in COSMIC database

This stage will generally need to be run just once to generate the database.
Note that this is build specific (GRCh37 implemented), and the dbSNP and COSMIC databases have regular updates.

Note also that need to get credentials to download the COSMIC database.  Username and password
are stored in the file `COSMIC_credentials.dat` (run `2_download_cosmic.sh` to create sample
file.  This file is not tracked in the repository.)

Also need to create short test VCF for Strelka test data.  This done for development to create
data which is distributed with this project.  This step is not generally run again.
