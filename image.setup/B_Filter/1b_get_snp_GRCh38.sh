# Download b150 version of all variant VCF, and retain only headers and first 5 columns of VCF

# GRCh37 All 
SRC="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/00-All.vcf.gz"
OUT="human_9606_b150_GRCh38p7.All.5col.vcf.gz"

bash ./get_snp.sh $SRC $OUT
