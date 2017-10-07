# somaticwrapper version 1.0 ##

Detect somatic variants from tumor and normal exome data

### Song Cao ###

SomaticWrapper pipeline is a fully automated and modular software package designed for detection of somatic variants from tumor and normal exome data. It works on LSF job scheduler and can run multiple jobs in parallel. It was developed from GenomeVIP. Multiple standard variant calling tools are included in the pipeline such as varscan, strelka, mutect and pindel.

Pipeline version: 1

Usage: perl somatic_calling_v1.1.pl  --srg --step --sre --rdir --ref 

rdir = full path of the folder holding files for this sequence run (user must provide)

srg = bam having read group or not: 1, yes and 0, no (default 1)

sre = re-run: 1, yes and 0, no  (default 0)

step = run this pipeline step by step. (user must provide)

ref = the human reference: 

with chr: /gscmnt/gc3027/dinglab/medseq/fasta/GRCh37V1/GRCh37-lite-chr_with_chrM.fa

without chr: /gscmnt/gc3027/dinglab/medseq/fasta/GRCh37/GRCh37-lite.fa 

mmy: /gscmnt/gc2737/ding/Reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.20170530.fa 

  [0]  Run all steps
 
  [1]  Run streka
 
  [2]  Run Varscan

  [3]  Parse streka result

  [4]  Parse VarScan result

  [5]  Run Pindel

  [6]  Parse Pindel

  [7]  Run VEP annotation

  [8]  Run mutect 
  
  [9]  Merge vcf files  

  [10] generate maf file

  [11] Merge maf together 
 
