# from http://snpeff.sourceforge.net/SnpSift.html#concordance

if [ -z $1 ] || [ -z $2 ]; then
echo Usage:  $0 A.vcf B.vcf 
exit
fi

if [ ! -f $1 ]; then
echo $1 does not exist
exit
fi

if [ ! -f $2 ]; then
echo $2 does not exist
exit
fi

echo Comparing $1 and $2

JAR="/usr/local/snpEff/SnpSift.jar"
#JAR="/gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar"
java -jar $JAR concordance -v $1 $2
