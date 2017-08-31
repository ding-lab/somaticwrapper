#!/bin/bash

#cosmic.dbsnp.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar
#cosmic.dbsnp.db = /gscuser/scao/gc3027/cosmic/CosmicCodingMuts.vcf
#cosmic.dbsnp.rawvcf = /gscmnt/gc2525/dinglab/rmashl/Software/bin/dbSNP/NCBI/snp142/GRCh37/00-All.brief.vcf
#cosmic.dbsnp.mode = filter
#cosmic.dbsnp.passfile  = /gscuser/scao/gc3027/cosmic/00-All.brief.pass.cosmic.vcf
#cosmic.dbsnp.dbsnpfile = /gscuser/scao/gc3027/cosmic/00-All.brief.not.pass.cosmic.vcf

DATD="/data/B_Filter"
JAR="/usr/local/snpEff/SnpSift.jar"

RAWVCF="$DATD/short.5col.vcf.bgz"
COSMICDB="$DATD/CosmicCodingMuts.vcf.bgz"



anno="$DATD/dbsnp_anno.vcf";

java -jar $JAR annotate -id $COSMICDB $RAWVCF > $anno

echo Written to $anno
exit

#my $cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} annotate -id $paras{'db'} $paras{'rawvcf'} > $anno";
#print "$cmd\n";
#system($cmd);
#checksize($anno, $paras{'rawvcf'});
#if( exists $paras{'mode'}  &&  $paras{'mode'} eq "filter" )  {
#$cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} filter -n \" (exists ID) & (ID =~ 'COSM' ) \" -f $anno > $paras{'passfile'}";
