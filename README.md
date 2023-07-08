
# somaticwrapper: tumor only pipeline, compute1 version 1.0  #

Detect somatic variants from tumor only data based on HG38 reference. SomaticWrapper pipeline is a fully automated and modular software package designed for detection of somatic variants from tumor exome data. It works on LSF job scheduler and can run multiple jobs in parallel. 

Use VEP85 for annotation 
## Install the third-party software ##


## Usage ##

Step1: Enter the directory where you downloaded somaticwrapper pipeline 

Step2: Type the coommand line: perl somaticwrapper.pl --rdir --ref --log --q --groupname --users --step 

rdir = full path of the folder holding files for this sequence run (user must provide)

log = full path of the folder for saving log file; usually upper folder of rdir

step run this pipeline step by step. (user must provide)

ref: the human reference: 

q: which queue for submitting job; research-hpc, ding-lab, long (default)

groupname = job group name: Format; users/groupname

users = user name for job group: Format; users/groupname

[0]  generate bams if input files are fastqs

[1]  Run Mutect2

[2]  Run filter Mutect2 result

[3]  Run parse Mutect2 result

[4] Generate maf file

[5] Generate merged maf file

## Example ##

perl somaticwrapper.pl --rdir /storage1/fs1/dinglab/Active/Projects/scao/gbm/GSAM --log /storage1/fs1/dinglab/Active/Projects/scao/gbm/GSAM.log --ref  /storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa --q general --users songcao --groupname SomaticWXS --step 2

## Contact ##

Song Cao, scao@wustl.edu 
