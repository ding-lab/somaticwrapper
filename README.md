
# Somaticwrapper version 2.2, compute1 #

Detect somatic variants from tumor and normal WGS/WXS data (HG38 reference). SomaticWrapper pipeline is a fully automated and modular software package designed for detection of somatic variants from tumor and normal exome data. It works on LSF job scheduler and can run multiple jobs in parallel. Multiple standard variant calling tools are included in the pipeline such as varscan2, strelka2, mutect1 and pindel. The final called variant can be found from dnp.annotated.maf for all variants and dnp.annotated.coding.maf for coding variants.  

SNV calls are intersecting results from 2 over 3 callers (Strelka2, Mutect1, and VarScan2).

Indel calls are called by 2 over 3 callers (Strelka2, Varscan2 and pindel). 

Improvements compared to version v2.1:

1. Remove indels > 100 nt before annotation

2. Fix false alarm for step 7

## Usage ##

Step 0: set environment for LSF job on compute1 by adding the following to ~/.bashrc file: 

export PATH=/storage1/fs1/songcao/Active/Software/anaconda3/bin:$PATH

export STORAGE2=/storage1/fs1/dinglab/Active
export SCRATCH2=/storage1/fs1/dinglab/

export STORAGE1=/storage1/fs1/songcao/Active
export SCRATCH1=/storage1/fs1/songcao/

export LSF_DOCKER_VOLUMES="$STORAGE1:$STORAGE1 $STORAGE2:$STORAGE2"

Then, start interactive queue environment for running jobs: 
bsub -G compute-dinglab -q dinglab-interactive -Is -a 'docker(scao/dailybox)' /bin/bash

Step1: Enter the directory where you downloaded somaticwrapper pipeline 

Step2: Type the coommand line: perl somaticwrapper.pl  --srg --sre --wgs --rdir --ref --log --q --mincovt --mincovn --minvaf --maxindsize --exonic --smg --groupname --users --step

rdir = full path of the folder holding files for this sequence run (user must provide)

log = full path of the folder for saving log file; usually upper folder of rdir

srg = bam having read group or not: 1, yes and 0, no (default 1)

sre = re-run: 1, yes and 0, no  (default 0)

wgs = 1 if it is wgs data and otherwise it is 0; If you want to output the maf for all variants, set exonic to 0

groupname = job group name: Format; users/groupname

users = user name for job group: Format; users/groupname

step run this pipeline step by step. (user must provide)

ref: the human reference: 

q: which queue for submitting job; research-hpc, ding-lab, long (default)

mincovt: minimum coverage for tumor: default >=14

mincovn: minimum coverage for normal: default >=8

minvaf: minimum somatic vaf: default >=0.05

maxindsize: default <= 100

exonic: output exonic region: 1 Yes, 0 No, Default Yes

smg: smg gene list that escapes the 0.05 vaf cut-off



[1]  Run streka 

[2]  Run Varscan 

[3]  Run Pindel 

[4]  Run mutect 

[5]  Parse mutect result 

[6]  Parse streka result 

[7]  Parse VarScan result 

[8]  Parse Pindel 

[9] QC vcf files

[10]  Merge vcf files   

[11] Generate maf file  

[12] Generate merged maf file 

[13] DNP annotation

[14] Clean unnecessary intermediate files

## Contact ##

Song Cao, scao@wustl.edu 
