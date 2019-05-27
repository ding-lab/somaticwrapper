# somaticwrapper version 1.4 ##

Detect somatic variants from tumor and normal WXS based on HG38 reference

### Song Cao ###

*you must enter the directory with somaticwrapper pipeline to submit the jobs

SomaticWrapper pipeline is a fully automated and modular software package designed for detection of somatic variants from tumor and normal exome data. It works on LSF job scheduler and can run multiple jobs in parallel. Multiple standard variant calling tools are included in the pipeline such as varscan, strelka, mutect and pindel.

Pipeline version: 1.4

Usage: perl somaticwrapper.pl  --srg --step --sre --rdir --ref --log --q --mincovt --mincovn --minvaf --maxindsize

rdir = full path of the folder holding files for this sequence run (user must provide)

log = full path of the folder for saving log file; usually upper folder of rdir

srg = bam having read group or not: 1, yes and 0, no (default 1)

sre = re-run: 1, yes and 0, no  (default 0)

step run this pipeline step by step. (user must provide)

ref: the human reference: 

q: which queue for submitting job; research-hpc, ding-lab, long (default)

mincovt: minimum coverage for tumor: default >=14

mincovn: minimum coverage for normal: default >=8

minvaf: minimum somatic vaf: default >=0.05

maxindsize: default <=100

hg38: /gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data/image.data/A_Reference/GRCh38.d1.vd1.fa

[0]  Run all steps 

[1]  Run streka 

[2]  Run Varscan 

[3]  Run Pindel 

[4]  Run mutect 

[5]  Parse mutect result 

[6]  Parse streka result 

[7]  Parse VarScan result 

[8]  Parse Pindel 

[9]  Merge vcf files   

[10] Generate maf file  

[11] Generate merged maf file 

[12] DNP annotation 
