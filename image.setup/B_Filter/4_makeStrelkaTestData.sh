# Create short test VCF for Strelka test data
# This is generated just for development purposes, and the resulting data
# is distributed in this project
# In most cases there is no need to run this step

# Strelka test data has reads in region 20:900-4000
# There are no snps in dbsnp VCF, but there is a large number of snps in region 20:1000900-1004000,
# so will create test dbSnP database by subtracting 1000000 from each position 

DATD="/data/image.data/B_Filter"
OUTD="../S_StrelkaDemoSetup/data"
mkdir -p $OUTD

DBSNP="human_9606_b150_GRCh37p13.All.5col.vcf.gz"
COSMIC="CosmicCodingMuts.grch37.v82.vcf.gz"

OFFSET=1000000
DBSNP_S="$OUTD/dbSnP.vcf.gz"
tabix $DATD/$DBSNP 20:1000900-1004000 | awk -v offset=$OFFSET 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2-offset,$3,$4,$5,$6,$7,$8}' | bgzip > $DBSNP_S
tabix -p vcf $DBSNP_S
echo Written to $DBSNP_S

#For Cosmic, the above region does not have anything, but the region below has a cluster of variants
#tabix CosmicCodingMuts.grch37.v82.vcf.gz 20:1099445-1100000 
#(subtract 1098545 from above to map to 900-4000)

OFFSET=1098545
COSMIC_S="$OUTD/COSMIC.vcf.gz"
tabix $DATD/$COSMIC 20:1099445-1100000 | awk -v offset=$OFFSET 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2-offset,$3,$4,$5,$6,$7,$8}' | bgzip > $COSMIC_S
tabix -p vcf $COSMIC_S
echo Written to $COSMIC_S

# TODO: confirm script below works

OUT="$OUTD/dbsnp-StrelkaDemo.noCOSMIC.vcf.gz"
bash ./make_variant_filter.sh $DBSNP_S $COSMIC_S $OUT


