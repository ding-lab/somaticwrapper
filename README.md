
# somaticwrapper: tumor only pipeline  #

Detect somatic variants from tumor only data based on HG38 reference. SomaticWrapper pipeline is a fully automated and modular software package designed for detection of somatic variants from tumor exome data. It works on LSF job scheduler and can run multiple jobs in parallel. 

## Install the third-party software ##


## Usage ##

Step1: Enter the directory where you downloaded somaticwrapper pipeline 

Step2: Type the coommand line: perl somaticwrapper.pl  --step --rdir --ref --log --q

rdir = full path of the folder holding files for this sequence run (user must provide)

log = full path of the folder for saving log file; usually upper folder of rdir

step run this pipeline step by step. (user must provide)

ref: the human reference: 

q: which queue for submitting job; research-hpc, ding-lab, long (default)

hg38: /gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data/image.data/A_Reference/GRCh38.d1.vd1.fa

[0]  Run fastq trimming

[1]  Generate bam

[2]  Run Varscan 


## Contact ##

Song Cao, scao@wustl.edu 
