
SAMPLE_NAME=$1

SAMPLE_NAME="TCGA-HI-TEST.GRCh37"

DATD="/data/data/$SAMPLE_NAME"
VEP="$DATD/merged/merged.VEP.vcf.vep"

CANCER_GENES="Cancer.Gene.List.dat"

# Write basic implementation of filter_vep analysis on test dataset
VEP="/data/data/$SAMPLE_NAME/merged/merged.VEP.vcf.vep"
echo Evaluating $VEP


while read GENE; do
# Skip comments 
[[ $GENE = \#* ]] && continue
#echo Evaluating $GENE in $VEP

COUNT=$(bash ./get_gene_count.sh $VEP $GENE)

if [ -z $COUNT ]; then
COUNT="0"
fi

BUILD=$(echo $SAMPLE_NAME | cut -f 2 -d .)
PARTICIPANT=$(echo $SAMPLE_NAME | cut -f 1 -d .)

printf "$SAMPLE_NAME\t$PARTICIPANT\t$BUILD\t$GENE\t$COUNT\n"

done <$CANCER_GENES

