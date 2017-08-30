#!/bin/bash

RUNDIR=/data/data/SWtest
VARSCAN_JAR=/usr/local/VarScan.v2.3.8.jar
SAMTOOLS=/usr/local/bin/samtools
#export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64

# these were exported.  Not clear that's necessary
JAVA_OPTS="-Xms256m -Xmx512m"
PATH=${JAVA_HOME}/bin:${PATH}

export LD_LIBRARY_PATH=${JAVA_HOME}/lib:${LD_LIBRARY_PATH}

echo Running in ${RUNDIR}/varscan
cd ${RUNDIR}/varscan

TMPBASE=./varscan.out.som
LOG=${TMPBASE}.log
snvoutbase=${TMPBASE}_snv
indeloutbase=${TMPBASE}_indel

BAMLIST=${RUNDIR}/varscan/bamfilelist.inp
if [ ! -e ${BAMLIST} ]
then
rm ${BAMLIST}
fi
echo "/data/data/SWtest/SWtest.N.bam" > ${BAMLIST}
echo "/data/data/SWtest/SWtest.T.bam" >> ${BAMLIST}

# with BAMLIST defined as above, ncols will always be 9
# getting rid of line below obviates dependency on bc package
## ncols=$(echo "3*( $(wc -l < $BAMLIST) +1)"|bc)
## echo ncols is $ncols   # its always 9

echo Writing to $LOG

REF="/data/A_Reference/demo20.fa"

SAMTOOLS_CMD="$SAMTOOLS mpileup -q 1 -Q 13 -B -f $REF -b ${BAMLIST} "

ARGS=" --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal 30 --min-coverage-tumor 22 --min-var-freq 0.08 --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp ${snvoutbase} --output-indel ${indeloutbase} "

JAVA_CMD="java ${JAVA_OPTS} -jar $VARSCAN_JAR somatic - ${TMPBASE} $ARGS"

# Note, we're removing the awk pipe below.  Does not make a difference in test case
# seems to be a filter for only lines with ncols columns?
#$SAMTOOLS_CMD | awk -v ncols=$ncols 'NF==ncols' | $JAVA_CMD &> ${LOG}

$SAMTOOLS_CMD | $JAVA_CMD &> ${LOG}

#Original command
#${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f  -b ${BAMLIST} | awk -v ncols=$ncols 'NF==ncols' | java ${JAVA_OPTS} -jar ${VARSCAN_DIR}/VarScan.jar somatic - ${TMPBASE} --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal 30 --min-coverage-tumor 22 --min-var-freq 0.08 --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp ${snvoutbase} --output-indel ${indeloutbase} &> ${LOG}
