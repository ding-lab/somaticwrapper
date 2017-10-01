# Usage: 
# get_gene_count.sh CCDC132
# This is equivalent to passing --filter "SYMBOL is CCDC132" to filter_vep

# NOTE: important to put quotes around input string on command line

VEP=$1
GENE=$2


# Based on:
# perl /usr/local/ensembl-vep/vep  --buffer_size 10000 --offline --cache --dir /data/D_VEP --assembly GRCh37 --fork 4 --format vcf --symbol -i /tmp/MHq97ukYyj -o /tmp/fYO5Wh2_Vu --force_overwrite  --fasta /data/A_Reference/Homo_sapiens_assembly19.fasta

#  This is a ad hoc parser
# grep -v "^##" /data/data/TCGA-A2-A0CM.GRCh37/merged/merged.VEP.vcf.vep | cut -f 14 | tr ';' '\n' | grep SYMBOL= | cut -f 2 -d = | sort | uniq -c | sort -nr | head

# Background: http://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html

FILTER_VEP="perl /usr/local/ensembl-vep/filter_vep"

#perl /usr/local/ensembl-vep/filter_vep -i /data/data/TCGA-HI-TEST/merged/merged.VEP.vcf.vep -filter "SYMBOL is CCDC132"

FILTER="SYMBOL is $GENE"

# By default writes to STDOUT
$FILTER_VEP -i $VEP -filter "$FILTER" -c

