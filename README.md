# somaticwrapper version 2.1 ##
### for HG38 ##

Detect coding somatic variants from tumor and normal WXS data

### Song Cao ###

### ********you must enter the directory with somaticwrapper pipeline to submit the jobs******* ###

SomaticWrapper pipeline is a fully automated and modular software package designed for detection of somatic variants from tumor and normal exome data. It works on LSF job scheduler and can run multiple jobs in parallel. Multiple standard variant calling tools are included in the pipeline such as varscan, strelka and pindel.

Pipeline version: 1.1

Usage: perl somatic_calling_v1.1.pl --srg --step --sre --rdir --ref --log --q --wgs --indsize 

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<srg> = bam having read group or not: 1, yes and 0, no (default 1)
<sre> = re-run: 1, yes and 0, no  (default 0)
<refname> = GRCh37 or Hg19
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
<q> which queue for submitting job; research-hpc, ding-lab, long (default)
<wgs> ==  1 for yes and 0 for no 
<indsize> = indel size < indsize; default indsize = 100
HG38: 



  [0]  Run all steps
 
  [1]  Run mutect
 
  [2]  Run Varscan

  [3] Run Pindel

  [4]  Parse mutect result

  [5]  Parse VarScan result

  [6]  Parse Pindel
  
  [7]  Merge vcf files  

  [8] generate maf file

  [9] Merge maf together 
 
