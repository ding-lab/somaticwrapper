# somaticwrapper version 1.1 ##

Detect somatic variants from tumor and normal WXS

### Song Cao ###

SomaticWrapper pipeline is a fully automated and modular software package designed for detection of somatic variants from tumor and normal exome data. It works on LSF job scheduler and can run multiple jobs in parallel. Multiple standard variant calling tools are included in the pipeline such as varscan, strelka, mutect and pindel.

Pipeline version: 1.1

$yellow     Usage: perl $0  --srg --step --sre --rdir --ref  --log --q 

$normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<srg> = bam having read group or not: 1, yes and 0, no (default 1)
<sre> = re-run: 1, yes and 0, no  (default 0)
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
<q> which queue for submitting job; research-hpc, ding-lab, long (default)

hg38: /gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data/image.data/A_Reference/GRCh38.d1.vd1.fa

$red         [0]  Run all steps
$green       [1]  Run streka
$green       [2]  Run Varscan
$green       [3]  Run Pindel
$green       [4]  Run mutect
$yellow      [5]  Parse mutect result
$yellow      [6]  Parse streka result
$yellow      [7]  Parse VarScan result
$yellow      [8]  Parse Pindel
$cyan        [9]  Merge vcf files  
$cyan        [10] Generate maf file 
$cyan        [11] Generate merged maf file
$normal
 
