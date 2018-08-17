##########################################################################################
### Author: Song Cao
### Author: Sunantha Sethuraman
### Description: Original somatic wrapper developed by Song adapted for Katmai by Sunantha
##########################################################################################

##########################################################################################
### Initialize perl
##########################################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $version = 1.0;
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";

##########################################################################################
### Usage
##########################################################################################

(my $usage = <<OUT) =~ s/\t+//g;

Somatic Wrapper: A variant calling pipeline 
Pipeline version: $version

This wrapper script can only handle hg19/GRCh37 build. It is NOT suitable for use with hg38/GRCh38.

$yellow     Usage: perl $0  --srg --step --sre --rdir --ref --refname --log --wgs 

$normal

<rdir> = Full path to sample location ### See directory structure below 
<log> = Full path to desired log directory ### See recommendation below
<srg> = BAM contains read groups: 1 for yes and 0 for no (default 1)
<sre> = Re-run: 1 for yes and 0 for no  (default 0)
<refname> = GRCh37 or Hg19 ## Case sensitive
<step> = Step number to run ## See steps below
<ref> == Full path to the reference fasta file # MUST be identical to the fasta used for alignment of BAMs, MUST contain an index file in the same directory
<wgs> ==  Whole genome sequencing: 1 for yes and 0 for no 

$normal

Directory structure: The samples should be stored in a certain directory format to ensure that Somatic Wrapper works as intended. 
If you have two samples with names "mysample1" and "exsample", the recommended directory format would be:

$cyan 

---desired location on your computer
	---data
		----mysample
			----mysample.T.bam
			----mysample.T.bam.bai
			----mysample.N.bam
			----mysample.N.bam.bai
		----exsample
			---exsample.T.bam
			---exsample.T.bam.bai
			---exsample.N.bam
			---exsample.N.bam.bai
	---log

Please note that the bam files MUST be indexed, with the index files (.bai files) named as per the convention shown above.
In the above examples, the "rdir" would be ".../desired location on your computer/data" and "log" would be ".../desired location on your computer/log".
Somatic Wrapper will attempt to consider as input ALL the sub-directories of "rdir", so please do not have any other subdirectories or files present inside "rdir".
All shell scripts generated in the process will be written to the log directory.


$red	[0] Run all steps  ### Does not work on Katmai at the moment
$green	[1] Run strelka
$green	[2] Run VarScan
$green	[3] Run Pindel
$yellow	[4] Run Parse Strelka
$yellow	[5] Run Parse VarScan
$yellow	[6] Run Parse Pindel
$red	[7] Run Merge VCF 
$red	[8] Run vcf2maf
$red	[9] Generate final report

$normal

OUT

##########################################################################################
### Set defaults and parse arguments
##########################################################################################

#__HELP (BOOLEAN, DEFAULTS TO NO-HELP)
my $help = 0;
#__DEFAULT NUMBER OF BINS IE (MUST BE INTEGER)
my $step_number = -1;
my $status_rg = 1;
my $status_rerun=0;
my $s_wgs=0;
#__FILE NAME (STRING, NO DEFAULT)
my $run_dir="";
my $log_dir="";
my $h37_REF="";
my $ref_name="";

#__PARSE COMMAND LINE
my $status = &GetOptions (
	"step=i" => \$step_number,
	"srg=i" => \$status_rg,
	"sre=i" => \$status_rerun,
	"rdir=s" => \$run_dir,
	"ref=s"  => \$h37_REF,
	"log=s"  => \$log_dir,
	"refname=s"  => \$ref_name,
	"wgs=i"  => \$s_wgs,
	"help" => \$help, 
	);
 
##########################################################################################
### Starting checks
##########################################################################################

### Check if all mandatory arguments are provided, else print help message

if ($help || $run_dir eq "" ||  $ref_name eq "" || $log_dir eq "" || $step_number<0) {
	print $usage;
    exit;
}

### Check if user provided step number is within the acceptable range, else print help message

die $usage unless ($step_number >=0)&&(($step_number <= 9));

### Check if run directory exists and if accessible, read sample list

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

##########################################################################################
### Assign scripts and log directory
##########################################################################################

### Assign a directory to write shell scripts to

if (! -d $log_dir."/somaticwrapper_scripts") {
	`mkdir $log_dir"/somaticwrapper_scripts"`;
}
my $job_files_dir = $log_dir."/somaticwrapper_scripts";

### Assign a directory to write log and stdout to

if (! -d $log_dir."/somaticwrapper_log") {
    `mkdir $log_dir"/somaticwrapper_log"`;
}
my $log_file_dir = $log_dir."/somaticwrapper_log";


##########################################################################################
### Set a few global variables
##########################################################################################

### Set script directory (this is the somatic wrapper path, not the path to shell scripts)

my $run_script_path = dirname(abs_path($0));
my $run_script_perl = "/usr/bin/perl ".$run_script_path."/";

### Set HOME and initialize a few empty variables

my $HOME = $ENV{HOME};
my $current_job_file = "";
my $sample_full_path = "";
my $sample_name = "";
my $sub_com = "";

### Print user provided arguments for each run

print "run dir=",$run_dir,"\n";
print "log dir=",$log_dir,"\n";
print "step num=",$step_number,"\n";
print "status rerun=",$status_rerun,"\n";
print "status readgroup=",$status_rg,"\n";


##########################################################################################
### Specify paths to programs on Katmai (hardcoded)
##########################################################################################


my $STRELKA_DIR="/diskmnt/Software/strelka_workflow-1.0.14/bin/bin";
my $strelka_confdir="/diskmnt/Software/Tools_for_SomaticWrapper/config_strelka";
my $VARSCAN_DIR="/diskmnt/Software/varscan-2.3.8";
my $varscan_confdir="/diskmnt/Software/Tools_for_SomaticWrapper/config_varscan";
my $f_exac="/diskmnt/Datasets/ExAC/r0.3.1/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
my $f_ref_annot="/diskmnt/Software/VEP_v85/cache/homo_sapiens/85_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
my $pindel="/diskmnt/Software/pindel-0.2.5b9-commit-b706fba/pindel/pindel";
my $PINDEL_DIR="/diskmnt/Software/pindel-0.2.5b9-commit-b706fba/pindel/";
my $SAMTOOLS_DIR="/diskmnt/Software/samtools-1.2/bin";
my $snpsift="/diskmnt/Software/snpEff_20150522/SnpSift.jar";
my $cosmicvcf="/diskmnt/Software/Tools_for_SomaticWrapper/00-All.brief.pass.cosmic.vcf";
my $bamreadcount="/diskmnt/Software/bam-readcount-0.7.4/mybuild/bin/bam-readcount";
my $vepscript="/diskmnt/Software/VEP_v85/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl";
my $vepcache="/diskmnt/Software/VEP_v85/cache";
my $picardexe="/diskmnt/Software/picard.jar";
my $gatk="/diskmnt/Software/GenomeAnalysisTK-3.7-0-gcfedb67.jar";
my $java_dir="/diskmnt/Software/jre1.8.0_121";
my $f_centromere="/diskmnt/Datasets/Centromeres/gc3015/dinglab/medseq/Jiayin_Germline_Project/PCGP/data/pindel-centromere-exclude.bed";


##########################################################################################
### Calling programs step-by-step
##########################################################################################

for (my $i=0;$i<@sample_dir_list;$i++) { 
	$sample_name = $sample_dir_list[$i];
	
	if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
		$sample_full_path = $run_dir."/".$sample_name;
		
		if (-d $sample_full_path) {
			print $yellow, "\nSubmitting job for the sample ",$sample_name, " ",$normal, "\n";
			$current_job_file="";
			
			if($step_number==0) {
				#&sub_strelka();
				#&sub_varscan();
				#&sub_pindel();
				#&sub_parse_strelka();
				#&sub_parse_varscan();
				#&sub_parse_pindel();
				#&sub_merge_vcf();
				#&sub_vcf_2_maf();
			}
			
			elsif ($step_number == 1) {
				&sub_strelka();
			}
			
			elsif ($step_number == 2) {
				&sub_varscan(1);
			}
				
			elsif ($step_number == 3) {
				&sub_pindel(1);
			}
				
			elsif ($step_number == 4) {
				&sub_parse_strelka(1);
			}
			
			elsif ($step_number == 5) {
				&sub_parse_varscan(1);
			}
			
			elsif ($step_number == 6) {
				&sub_parse_pindel(1);
			}
			
			elsif ($step_number == 7) {
				&sub_merge_vcf(1);
			}
			
			elsif ($step_number == 8) {
				&sub_vcf_2_maf(1);
			}
			
			elsif ($step_number == 9) {
				&sub_gen_report(1);
			}
		}
	}
}

##########################################################################################
### Step 1: Run Strelka
##########################################################################################

sub sub_strelka{
	
	$current_job_file = "j1_strelka_".$sample_name.".sh"; 
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
	
	### Check if input BAM exist
	
	if (! -e $IN_bam_T) {
		print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_T does not exist", $normal, "\n\n";
	}
	
	if (! -e $IN_bam_N) {
		print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		 die "Died because $IN_bam_N does not exist", $normal, "\n\n";
	}

	### Check that the input BAM are not empty
	
	if (! -s $IN_bam_T) {
		print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_T is empty!", $normal, "\n\n";
	}
	
	if (! -s $IN_bam_N) {#make sure input fasta file is not empty
		print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_N is empty!", $normal, "\n\n";
	}
	
	### Specify error and log files
	
	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;
	
	### Print out the shell script (job script)

	open(STRELKA, ">$job_files_dir/$current_job_file") or die $!;
	
	print STRELKA "#!/bin/bash\n";
	print STRELKA "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
	print STRELKA "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
	print STRELKA "myRUNDIR=".$sample_full_path."/strelka\n";
	print STRELKA "STATUSDIR=".$sample_full_path."/status\n";
	print STRELKA "RESULTSDIR=".$sample_full_path."/results\n";
	print STRELKA "SG_DIR=".$sample_full_path."/strelka\n"; 
	print STRELKA "RUNDIR=".$sample_full_path."\n";
	print STRELKA "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
	print STRELKA "STRELKA_VCF=".$sample_full_path."/strelka/strelka_out/results/passed.somatic.snvs.vcf"."\n";   
	print STRELKA "CONFDIR=".$strelka_confdir."\n";
	print STRELKA "TASK_STATUS=".$sample_full_path."/strelka/strelka_out/task.complete"."\n";
	print STRELKA "export SAMTOOLS_DIR=".$SAMTOOLS_DIR."\n";
	print STRELKA "export JAVA_HOME=$java_dir\n";
	print STRELKA "export JAVA_OPTS=\"-Xmx10g\"\n";
	print STRELKA "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print STRELKA "if [ ! -d \${myRUNDIR} ]\n";
	print STRELKA "then\n";
	print STRELKA "mkdir \${myRUNDIR}\n";
	print STRELKA "fi\n";
	print STRELKA "if [ $status_rerun -eq 1 ]\n";
	print STRELKA "  then\n";
	print STRELKA "rm \${TASK_STATUS}\n";
	print STRELKA "fi\n";
	print STRELKA "if [ ! -f \${STRELKA_VCF} ]\n";
	print STRELKA "  then\n";
	print STRELKA "rm \${TASK_STATUS}\n";
	print STRELKA "fi\n";
	print STRELKA "if [ ! -f  \${TASK_STATUS} ]\n";
	print STRELKA "then\n";
	print STRELKA "if [ -d \${STRELKA_OUT} ]\n";
	print STRELKA "then\n";
	print STRELKA "rm -rf \${STRELKA_OUT}\n";
	print STRELKA "fi\n";
	print STRELKA "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n"; 
	print STRELKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
	print STRELKA "else\n";
	print STRELKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
	print STRELKA "fi\n";
	print STRELKA ". $run_script_path/set_envvars\n";
	print STRELKA "if [ $s_wgs -eq 1 ]\n";
	print STRELKA "  then\n";
	print STRELKA "   ".$STRELKA_DIR."/configureStrelkaWorkflow.pl --normal \$NBAM --tumor \$TBAM --ref ". $h37_REF." --config $run_script_path/strelka.ini.wgs --output-dir \$STRELKA_OUT\n";
	print STRELKA "else\n";	
	print STRELKA "   ".$STRELKA_DIR."/configureStrelkaWorkflow.pl --normal \$NBAM --tumor \$TBAM --ref ". $h37_REF." --config $run_script_path/strelka.ini --output-dir \$STRELKA_OUT\n";
	print STRELKA "fi\n";
	print STRELKA "cd \$STRELKA_OUT\n";
	print STRELKA "make -j 16\n";
	print STRELKA "touch \${TASK_STATUS}\n";
	print STRELKA "fi\n";
	close STRELKA;
	
	my $sh_file=$job_files_dir."/".$current_job_file;
	
	### Submit job to background

	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}


##########################################################################################
### Step 2: Run VarScan
##########################################################################################

sub sub_varscan{

	$current_job_file = "j2_varscan_".$sample_name.".sh";
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
	
	### Check if input BAM exist
	
	if (! -e $IN_bam_T) {
		print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_T does not exist", $normal, "\n\n";
	}
	
	if (! -e $IN_bam_N) {
		 print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_N does not exist", $normal, "\n\n";
	}

	### Check that the input BAM are not empty
	
	if (! -s $IN_bam_T) {
		print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_T is empty!", $normal, "\n\n";
	}
	
	if (! -s $IN_bam_N) {#make sure input fasta file is not empty
		print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_N is empty!", $normal, "\n\n";
	}
	
	### Specify error and log files
	
	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;
	
	### Print out the shell script (job script)
	
	open(VARSCAN, ">$job_files_dir/$current_job_file") or die $!;
	
	print VARSCAN "#!/bin/bash\n";
	print VARSCAN "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
	print VARSCAN "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
	print VARSCAN "myRUNDIR=".$sample_full_path."/varscan\n";
	print VARSCAN "STATUSDIR=".$sample_full_path."/status\n";
	print VARSCAN "RESULTSDIR=".$sample_full_path."/varscan_results\n";
	print VARSCAN "RUNDIR=".$sample_full_path."\n";
	print VARSCAN "CONFDIR=".$varscan_confdir."\n";
	print VARSCAN "export VARSCAN_DIR=".$VARSCAN_DIR."\n";
	print VARSCAN "export SAMTOOLS_DIR=".$SAMTOOLS_DIR."\n";
	print VARSCAN "export JAVA_HOME=$java_dir\n";
	print VARSCAN "export JAVA_OPTS=\"-Xmx10g\"\n";
	print VARSCAN "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print VARSCAN "if [ ! -d \${myRUNDIR} ]\n";
	print VARSCAN "then\n";
	print VARSCAN "mkdir \${myRUNDIR}\n";
	print VARSCAN "fi\n";
	print VARSCAN "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
	print VARSCAN "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
	print VARSCAN "else\n";
	print VARSCAN "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
	print VARSCAN "fi\n";
	print VARSCAN "put_cmd=\"ln -s\"\n";
	print VARSCAN "del_cmd=\"rm -f\"\n";
	print VARSCAN "del_local=\"rm -f\"\n";
	print VARSCAN "statfile=complete.vs_som_snvindels\n";
	print VARSCAN "localstatus=\${myRUNDIR}\/status\/\${statfile}\n";
	print VARSCAN "if [ ! -d \${myRUNDIR}\/status ]\n";
	print VARSCAN "then\n";
	print VARSCAN "mkdir \${myRUNDIR}\/status\n";
	print VARSCAN "fi\n";
	print VARSCAN "if [ $status_rerun -eq 1 ]\n";
	print VARSCAN "  then\n";
	print VARSCAN "rm \${localstatus}\n";
	print VARSCAN "fi\n";
	print VARSCAN "if [ ! -f  \${localstatus} ]\n";
	print VARSCAN "then\n";
	print VARSCAN "cd \${RUNDIR}/varscan\n";
	print VARSCAN "TMPBASE=.\/varscan.out.som\n";
	print VARSCAN "LOG=\${TMPBASE}.log\n";
	print VARSCAN "snvoutbase=\${TMPBASE}_snv\n";
	print VARSCAN "indeloutbase=\${TMPBASE}_indel\n";
	print VARSCAN "BAMLIST=\${RUNDIR}/varscan/bamfilelist.inp\n";
	print VARSCAN "if [ ! -e \${BAMLIST} ]\n";
	print VARSCAN "then\n";
	print VARSCAN "rm \${BAMLIST}\n";
	print VARSCAN "fi\n";
	print VARSCAN "echo \"$IN_bam_N\" > \${BAMLIST}\n"; 
	print VARSCAN "echo \"$IN_bam_T\" >> \${BAMLIST}\n";  
	print VARSCAN "ncols=9\n"; 
	print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somatic - \${TMPBASE} --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal 20 --min-coverage-tumor 20 --min-var-freq 0.05 --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp \${snvoutbase} --output-indel \${indeloutbase} &> \${LOG}\n";
	print VARSCAN '      grep "Error occurred during initialization of VM" ${LOG}',"\n";# one possible blast error (see the end of this script). 
	print VARSCAN '      CHECK=$?',"\n";
	print VARSCAN '      while [ ${CHECK} -eq 0 ]',"\n";
	print VARSCAN "      do\n";
	print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somatic - \${TMPBASE} --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal 20 --min-coverage-tumor 20 --min-var-freq 0.05 --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp \${snvoutbase} --output-indel \${indeloutbase} &> \${LOG}\n";
	print VARSCAN '      grep "Error occurred during initialization of VM" ${LOG}',"\n";
	print VARSCAN '          CHECK=$?',"\n";
	print VARSCAN "      done\n";
	print VARSCAN "touch \${localstatus}\n";
	print VARSCAN "fi\n";
	close VARSCAN;
	
	my $sh_file=$job_files_dir."/".$current_job_file;
	
	### Submit job to background
	
	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}

##########################################################################################
### Step 3: Run Pindel
##########################################################################################

sub sub_pindel{
	
	$current_job_file = "j3_pindel_".$sample_name.".sh";  
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
	
	### Check if input BAM exist
	
	if (! -e $IN_bam_T) {
		print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_T does not exist", $normal, "\n\n";
	}
	
	if (! -e $IN_bam_N) {
		print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_N does not exist", $normal, "\n\n";
	}

	### Check that the input BAM are not empty
	
	if (! -s $IN_bam_T) {
		print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_T is empty!", $normal, "\n\n";
	}
	
	if (! -s $IN_bam_N) {#make sure input fasta file is not empty
		print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		die "Died because $IN_bam_N is empty!", $normal, "\n\n";
	}
	
	### Specify error and log files
	
	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;
	
	### Print out the shell script (job script)
	 
	open(PINDEL, ">$job_files_dir/$current_job_file") or die $!;
	
	print PINDEL "#!/bin/bash\n";
	print PINDEL "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
	print PINDEL "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
	print PINDEL "myRUNDIR=".$sample_full_path."/pindel\n"; 
	print PINDEL "CONFIG=\${myRUNDIR}"."/".$sample_name.".config\n";
	print PINDEL "statfile=complete.pindel\n";
	print PINDEL "localstatus=\${myRUNDIR}\/status\/\${statfile}\n";
	print PINDEL "if [ ! -d \${myRUNDIR} ]\n";
	print PINDEL "then\n";
	print PINDEL "mkdir \${myRUNDIR}\n";
	print PINDEL "fi\n"; 
	print PINDEL "if [ ! -d \${myRUNDIR}\/status ]\n";
	print PINDEL "then\n";
	print PINDEL "mkdir \${myRUNDIR}\/status\n";
	print PINDEL "fi\n";
	print PINDEL "echo \"$IN_bam_T\t500\t$sample_name.T\" > \${CONFIG}\n";
	print PINDEL "echo \"$IN_bam_N\t500\t$sample_name.N\" >> \${CONFIG}\n";
	print PINDEL "for chr in {1..22} X \n";
	print PINDEL "do \n";
	print PINDEL "nohup $pindel -T 4 -c \$chr -f $h37_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"."_\${chr}"." -m 6 -w 1 -J $f_centromere > \${myRUNDIR}"."/$sample_name"."_\${chr}_pindel.log"." & \n" ;
	print PINDEL "done \n";
	close PINDEL;
	
	my $sh_file=$job_files_dir."/".$current_job_file;
	
	### Submit job to background

	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}

##########################################################################################
### Step 4: Run Parse Strelka
##########################################################################################

sub sub_parse_strelka{

	$current_job_file = "j4_parse_strelka_".$sample_name.".sh";
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

	### Specify error and log files

	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;

	### Print out the shell script (job script)

	open(STRELKAP, ">$job_files_dir/$current_job_file") or die $!;

	print STRELKAP "#!/bin/bash\n";
	print STRELKAP "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
	print STRELKAP "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
	print STRELKAP "myRUNDIR=".$sample_full_path."/strelka\n";
	print STRELKAP "STATUSDIR=".$sample_full_path."/status\n";
	print STRELKAP "RESULTSDIR=".$sample_full_path."/results\n";
	print STRELKAP "SG_DIR=".$sample_full_path."/strelka\n";
	print STRELKAP "RUNDIR=".$sample_full_path."\n";
	print STRELKAP "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
	print STRELKAP "export SAMTOOLS_DIR=".$SAMTOOLS_DIR."\n";
	print STRELKAP "export VARSCAN_DIR=".$VARSCAN_DIR."\n";
	print STRELKAP "export JAVA_HOME=$java_dir\n";
	print STRELKAP "export JAVA_OPTS=\"-Xmx10g\"\n";
	print STRELKAP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print STRELKAP "cat > \${myRUNDIR}/strelka_out/results/strelka_dbsnp_filter.snv.input <<EOF\n";
	print STRELKAP "strelka.dbsnp.snv.annotator=".$snpsift."\n";
	print STRELKAP "strelka.dbsnp.snv.db=".$cosmicvcf."\n";
	print STRELKAP "strelka.dbsnp.snv.rawvcf = ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
	print STRELKAP "strelka.dbsnp.snv.mode = filter\n";
	print STRELKAP "strelka.dbsnp.snv.passfile  = ./strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
	print STRELKAP "strelka.dbsnp.snv.dbsnpfile = ./strelka.somatic.snv.all.gvip.dbsnp_present.vcf\n";
	print STRELKAP "EOF\n";
	print STRELKAP "cat > \${myRUNDIR}/strelka_out/results/strelka_dbsnp_filter.indel.input <<EOF\n";
	print STRELKAP "strelka.dbsnp.indel.annotator=".$snpsift."\n";
	print STRELKAP "strelka.dbsnp.indel.db = ".$cosmicvcf."\n";
	print STRELKAP "strelka.dbsnp.indel.rawvcf = ./strelka.somatic.indel.strlk_pass.gvip.vcf\n";
	print STRELKAP "strelka.dbsnp.indel.mode = filter\n";
	print STRELKAP "strelka.dbsnp.indel.passfile  = ./strelka.somatic.indel.all.gvip.dbsnp_pass.vcf\n";
	print STRELKAP "strelka.dbsnp.indel.dbsnpfile = ./strelka.somatic.indel.all.gvip.dbsnp_present.vcf\n";
	print STRELKAP "EOF\n";
	print STRELKAP "FP_BAM=\`awk \'{if(NR==1){print \$1}}\' \${RUNDIR}/varscan/bamfilelist.inp\`\n";
	print STRELKAP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_fpfilter.snv.input <<EOF\n";
	print STRELKAP "strelka.fpfilter.snv.bam_readcount =".$bamreadcount."\n";
	print STRELKAP "strelka.fpfilter.snv.bam_file = $IN_bam_T\n";
	print STRELKAP "strelka.fpfilter.snv.REF = $h37_REF\n"; 
	print STRELKAP "strelka.fpfilter.snv.variants_file = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
	print STRELKAP "strelka.fpfilter.snv.passfile = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_pass.vcf\n";
	print STRELKAP "strelka.fpfilter.snv.failfile = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_fail.vcf\n";
	print STRELKAP "strelka.fpfilter.snv.rc_in = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.in.vcf\n";
	print STRELKAP "strelka.fpfilter.snv.rc_out = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.out.vcf\n";
	print STRELKAP "strelka.fpfilter.snv.fp_out = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp.out.vcf\n";
	print STRELKAP "strelka.fpfilter.snv.min_mapping_qual = 0\n";
	print STRELKAP "strelka.fpfilter.snv.min_base_qual = 15\n";
	print STRELKAP "strelka.fpfilter.snv.min_num_var_supporting_reads = 4\n";
	print STRELKAP "strelka.fpfilter.snv.min_var_allele_freq = 0.05\n";
	print STRELKAP "strelka.fpfilter.snv.min_avg_rel_read_position = 0.10\n";
	print STRELKAP "strelka.fpfilter.snv.min_avg_rel_dist_to_3prime_end = 0.10\n";
	print STRELKAP "strelka.fpfilter.snv.min_var_strandedness = 0.01\n";
	print STRELKAP "strelka.fpfilter.snv.min_allele_depth_for_testing_strandedness = 5\n";
	print STRELKAP "strelka.fpfilter.snv.min_ref_allele_avg_base_qual = 30\n";
	print STRELKAP "strelka.fpfilter.snv.min_var_allele_avg_base_qual = 30\n";
	print STRELKAP "strelka.fpfilter.snv.max_rel_read_length_difference = 0.25\n";
	print STRELKAP "strelka.fpfilter.snv.max_mismatch_qual_sum_for_var_reads = 150\n";
	print STRELKAP "strelka.fpfilter.snv.max_avg_mismatch_qual_sum_difference = 150\n";
	print STRELKAP "strelka.fpfilter.snv.min_ref_allele_avg_mapping_qual = 30\n";
	print STRELKAP "strelka.fpfilter.snv.min_var_allele_avg_mapping_qual = 30\n";
	print STRELKAP "strelka.fpfilter.snv.max_avg_mapping_qual_difference = 50\n";
	print STRELKAP "EOF\n";
	print STRELKAP "cd \${STRELKA_OUT}/results\n";
	print STRELKAP "     ".$run_script_perl."genomevip_label.pl Strelka ./all.somatic.snvs.vcf ./strelka.somatic.snv.all.gvip.vcf\n";
	print STRELKAP "     ".$run_script_perl."genomevip_label.pl Strelka ./all.somatic.indels.vcf ./strelka.somatic.indel.all.gvip.vcf\n";
	print STRELKAP "     ".$run_script_perl."genomevip_label.pl Strelka ./passed.somatic.snvs.vcf ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
	print STRELKAP "     ".$run_script_perl."genomevip_label.pl Strelka ./passed.somatic.indels.vcf ./strelka.somatic.indel.strlk_pass.gvip.vcf\n";
	print STRELKAP "strelkasnvout=\${STRELKA_OUT}/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
	print STRELKAP "strelkaindelout=\${STRELKA_OUT}/results/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf\n";
	print STRELKAP "statfile=complete.strelka_parser\n";
	print STRELKAP "localstatus=\${STRELKA_OUT}\/\${statfile}\n";
	print STRELKAP "if [ $status_rerun -eq 1 ]\n";
	print STRELKAP "  then\n";
	print STRELKAP "rm \${localstatus}\n";
	print STRELKAP "fi\n";
	print STRELKAP "if [ ! -f  \${localstatus} ]\n";
	print STRELKAP "then\n"; 
	print STRELKAP "     ".$run_script_perl."dbsnp_filter.pl ./strelka_dbsnp_filter.snv.input\n";
	print STRELKAP "     ".$run_script_perl."dbsnp_filter.pl ./strelka_dbsnp_filter.indel.input\n";
	print STRELKAP "     ".$run_script_perl."snv_filter.pl ./strelka_fpfilter.snv.input\n";
	print STRELKAP '      grep "Error occurred during initialization of VM" ${strelkasnvout}',"\n";
	print STRELKAP '      CHECK=$?',"\n";
	print STRELKAP '      while [ ${CHECK} -eq 0 ]',"\n";
	print STRELKAP "      do\n";
	print STRELKAP "     ".$run_script_perl."dbsnp_filter.pl ./strelka_dbsnp_filter.snv.input\n";
	print STRELKAP "     ".$run_script_perl."dbsnp_filter.pl ./strelka_dbsnp_filter.indel.input\n";
	print STRELKAP "     ".$run_script_perl."snv_filter.pl ./strelka_fpfilter.snv.input\n";  
	print STRELKAP '      grep "Error occurred during initialization of VM" ${strelkasnvout}',"\n";
	print STRELKAP '          CHECK=$?',"\n";
	print STRELKAP "      done\n";
	print STRELKAP "touch \${localstatus}\n";
	print STRELKAP "  fi\n";
	close STRELKAP;

	my $sh_file=$job_files_dir."/".$current_job_file;

	### Submit job to background
	
	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}

##########################################################################################
### Step 5: Run Parse VarScan
##########################################################################################

sub sub_parse_varscan{

	$current_job_file = "j5_parse_varscan_".$sample_name.".sh";
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
	
	### Specify error and log files

	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;
	
	### Print out the shell script (job script)

	open(VARSCANP, ">$job_files_dir/$current_job_file") or die $!;

	print VARSCANP "#!/bin/bash\n";
	print VARSCANP "scr_t0=\`date \+\%s\`\n";
	print VARSCANP "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
	print VARSCANP "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
	print VARSCANP "myRUNDIR=".$sample_full_path."/varscan\n";
	print VARSCANP "STATUSDIR=".$sample_full_path."/status\n";
	print VARSCANP "RUNDIR=".$sample_full_path."\n";
	print VARSCANP "export VARSCAN_DIR=".$VARSCAN_DIR."\n";
	print VARSCANP "export SAMTOOLS_DIR=".$SAMTOOLS_DIR."\n";
	print VARSCANP "export JAVA_HOME=$java_dir\n";
	print VARSCANP "export JAVA_OPTS=\"-Xmx10g\"\n";
	print VARSCANP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_dbsnp_filter.snv.input <<EOF\n";
	print VARSCANP "varscan.dbsnp.snv.annotator =".$snpsift."\n";
	print VARSCANP "varscan.dbsnp.snv.db =".$cosmicvcf."\n";
	print VARSCANP "varscan.dbsnp.snv.rawvcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.snv.mode = filter\n";
	print VARSCANP "varscan.dbsnp.snv.passfile  = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.snv.dbsnpfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf\n";
	print VARSCANP "EOF\n";
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_dbsnp_filter.indel.input <<EOF\n";
	print VARSCANP "varscan.dbsnp.indel.annotator =".$snpsift."\n";
	print VARSCANP "varscan.dbsnp.indel.db =".$cosmicvcf."\n";
	print VARSCANP "varscan.dbsnp.indel.rawvcf = ./varscan.out.som_indel.gvip.Somatic.hc.vcf\n";
	print VARSCANP "varscan.dbsnp.indel.mode = filter\n";
	print VARSCANP "varscan.dbsnp.indel.passfile  = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.indel.dbsnpfile = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_present.vcf\n";
	print VARSCANP "EOF\n";
	print VARSCANP "FP_BAM=\`awk \'{if(NR==2){print \$1}}\' \${RUNDIR}/varscan/bamfilelist.inp\`\n";	
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_fpfilter.somatic.snv.input <<EOF\n";
	print VARSCANP "varscan.fpfilter.snv.bam_readcount =".$bamreadcount."\n";
	print VARSCANP "varscan.fpfilter.snv.bam_file = \${FP_BAM}\n";
	print VARSCANP "varscan.fpfilter.snv.REF = $h37_REF\n";
	print VARSCANP "varscan.fpfilter.snv.variants_file = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.fpfilter.snv.passfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.fp_pass.vcf\n";
	print VARSCANP "varscan.fpfilter.snv.failfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.fp_fail.vcf\n";
	print VARSCANP "varscan.fpfilter.snv.min_mapping_qual = 0\n";
	print VARSCANP "varscan.fpfilter.snv.min_base_qual = 15\n";
	print VARSCANP "varscan.fpfilter.snv.min_num_var_supporting_reads = 4\n";
	print VARSCANP "varscan.fpfilter.snv.min_var_allele_freq = 0.05\n";
	print VARSCANP "varscan.fpfilter.snv.min_avg_rel_read_position = 0.10\n";
	print VARSCANP "varscan.fpfilter.snv.min_avg_rel_dist_to_3prime_end = 0.10\n";
	print VARSCANP "varscan.fpfilter.snv.min_var_strandedness = 0.01\n";
	print VARSCANP "varscan.fpfilter.snv.min_allele_depth_for_testing_strandedness = 5\n";
	print VARSCANP "varscan.fpfilter.snv.min_ref_allele_avg_base_qual = 30\n";
	print VARSCANP "varscan.fpfilter.snv.min_var_allele_avg_base_qual = 30\n";
	print VARSCANP "varscan.fpfilter.snv.max_rel_read_length_difference = 0.25\n";
	print VARSCANP "varscan.fpfilter.snv.max_mismatch_qual_sum_for_var_reads = 150\n";
	print VARSCANP "varscan.fpfilter.snv.max_avg_mismatch_qual_sum_difference = 100\n";
	print VARSCANP "varscan.fpfilter.snv.min_ref_allele_avg_mapping_qual = 30\n";
	print VARSCANP "varscan.fpfilter.snv.min_var_allele_avg_mapping_qual = 30\n";
	print VARSCANP "varscan.fpfilter.snv.max_avg_mapping_qual_difference = 50\n";
	print VARSCANP "EOF\n";
	print VARSCANP "if [ ! -d \${myRUNDIR} ]\n";
	print VARSCANP "then\n";
	print VARSCANP "mkdir \${myRUNDIR}\n";
	print VARSCANP "fi\n";
	print VARSCANP "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
	print VARSCANP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
	print VARSCANP "else\n";
	print VARSCANP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
	print VARSCANP "fi\n";
	print VARSCANP "put_cmd=\"ln -s\"\n";
	print VARSCANP "del_cmd=\"rm -f\"\n";
	print VARSCANP "del_local=\"rm -f\"\n";
	print VARSCANP "statfile=complete.vs_som_parser\n";
	print VARSCANP "localstatus=\${RUNDIR}\/status\/\${statfile}\n";
	print VARSCANP "if [ ! -d \${myRUNDIR}\/status ]\n";
	print VARSCANP "then\n";
	print VARSCANP "mkdir \${myRUNDIR}\/status\n";
	print VARSCANP "fi\n";
	print VARSCANP "if [ $status_rerun -eq 1 ]\n";
	print VARSCANP "  then\n";
	print VARSCANP "rm \${localstatus}\n";
	print VARSCANP "fi\n";
	print VARSCANP "if [ ! -f  \${localstatus} ]\n";
	print VARSCANP "then\n";
	print VARSCANP "cd \${RUNDIR}/varscan\n";
	print VARSCANP "TMPBASE=.\/varscan.out.som\n";
	print VARSCANP "LOG=\${TMPBASE}.log\n";
	print VARSCANP "snvoutbase=\${TMPBASE}_snv\n";
	print VARSCANP "indeloutbase=\${TMPBASE}_indel\n";
	print VARSCANP "cd \${RUNDIR}/varscan\n";
	print VARSCANP "     ".$run_script_perl."genomevip_label.pl VarScan \${snvoutbase}.vcf  \${snvoutbase}.gvip.vcf\n";
	print VARSCANP "     ".$run_script_perl."genomevip_label.pl VarScan \${indeloutbase}.vcf \${indeloutbase}.gvip.vcf\n";
	print VARSCANP "echo \'APPLYING PROCESS FILTER TO SOMATIC SNVS:\' &>> \${LOG}\n";
	print VARSCANP "mysnvorig=./\${snvoutbase}.gvip.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar processSomatic \${mysnvorig} --min-tumor-freq 0.05 --max-normal-freq 0.05 --p-value 0.05  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_perl."extract_somatic_other.pl <  \${mysnvorig}  > \${mysnvorig/%vcf/other.vcf}\n";
	print VARSCANP "for kk in Somatic Germline LOH ; do\n";
	print VARSCANP "thisorig=\${mysnvorig/%vcf/\$kk.vcf}\n";
	print VARSCANP "thispass=\${mysnvorig/%vcf/\$kk.hc.vcf}\n";
	print VARSCANP "thisfail=\${mysnvorig/%vcf/\$kk.lc.vcf}\n";
	print VARSCANP "done\n";
	print VARSCANP "echo \'APPLYING PROCESS FILTER TO SOMATIC INDELS:\' &>> \$LOG\n";
	print VARSCANP "myindelorig=./\$indeloutbase.gvip.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar processSomatic \${myindelorig}   --min-tumor-freq  0.05   --max-normal-freq  0.05   --p-value  0.05  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_perl."extract_somatic_other.pl <  \${myindelorig}  >  \${myindelorig/%vcf/other.vcf}\n";
	print VARSCANP "for kk in Somatic Germline LOH ; do\n";
	print VARSCANP "thisorig=\${myindelorig/%vcf/\$kk.vcf}\n";
	print VARSCANP "thispass=\${myindelorig/%vcf/\$kk.hc.vcf}\n";
	print VARSCANP "thisfail=\${myindelorig/%vcf/\$kk.lc.vcf}\n";
	print VARSCANP "done\n";
	print VARSCANP "echo \'APPLYING SOMATIC FILTER:\' &>> \${LOG}\n";
	print VARSCANP "thissnvorig=\${snvoutbase}.gvip.Somatic.hc.vcf\n";
	print VARSCANP "myindelorig=\${indeloutbase}.gvip.vcf\n";
	print VARSCANP "thissnvpass=\${snvoutbase}.gvip.Somatic.hc.somfilter_pass.vcf\n";
	print VARSCANP "thissnvfail=\${snvoutbase}.gvip.Somatic.hc.somfilter_fail.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somaticFilter  ./\${thissnvorig} --min-coverage  20   --min-reads2  4   --min-strands2  1   --min-avg-qual  20   --min-var-freq  0.05 --p-value  0.05   --indel-file  ./\${myindelorig} --output-file  ./\${thissnvpass}  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_perl."dbsnp_filter.pl  \${RUNDIR}/varscan/vs_dbsnp_filter.snv.input\n";
	print VARSCANP "     ".$run_script_perl."dbsnp_filter.pl \${RUNDIR}/varscan/vs_dbsnp_filter.indel.input\n";
	print VARSCANP "touch \${localstatus}\n";
	print VARSCANP "fi\n";
	close VARSCANP; 

	my $sh_file=$job_files_dir."/".$current_job_file;
	
	### Submit job to background

	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}


##########################################################################################
### Step 6: Run Parse Pindel
##########################################################################################

sub sub_parse_pindel {

	$current_job_file = "j6_parse_pindel_".$sample_name.".sh";

	### Specify error and log files
	
	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;
	
	### Print out the shell script (job script)
	
	open(PINDELP, ">$job_files_dir/$current_job_file") or die $!;
	
	print PINDELP "#!/bin/bash\n";
	print PINDELP "RUNDIR=".$sample_full_path."\n";
	print PINDELP "cat > \${RUNDIR}/pindel/pindel_filter.input <<EOF\n";
	print PINDELP "pindel.filter.pindel2vcf = $PINDEL_DIR/pindel2vcf\n";
	print PINDELP "pindel.filter.variants_file = \${RUNDIR}/pindel/pindel.out.raw\n";
	print PINDELP "pindel.filter.REF = $h37_REF\n";
	print PINDELP "pindel.filter.date = 000000\n";
	print PINDELP "pindel.filter.heterozyg_min_var_allele_freq = 0.2\n";
	print PINDELP "pindel.filter.homozyg_min_var_allele_freq = 0.8\n";
	print PINDELP "pindel.filter.mode = somatic\n";
	print PINDELP "pindel.filter.apply_filter = true\n";
	print PINDELP "pindel.filter.somatic.min_coverages = 10\n";
	print PINDELP "pindel.filter.somatic.min_var_allele_freq = 0.05\n";
	print PINDELP "pindel.filter.somatic.require_balanced_reads = \"true\"\n";
	print PINDELP "pindel.filter.somatic.remove_complex_indels = \"true\"\n";
	print PINDELP "pindel.filter.somatic.max_num_homopolymer_repeat_units = 6\n";
	print PINDELP "EOF\n";
	print PINDELP "cat > \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input <<EOF\n";
	print PINDELP "pindel.dbsnp.indel.annotator =".$snpsift."\n";
	print PINDELP "pindel.dbsnp.indel.db =".$cosmicvcf."\n";
	print PINDELP "pindel.dbsnp.indel.rawvcf = ./pindel.out.current_final.gvip.Somatic.vcf\n";
	print PINDELP "pindel.dbsnp.indel.mode = filter\n";
	print PINDELP "pindel.dbsnp.indel.passfile  = ./pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
	print PINDELP "pindel.dbsnp.indel.dbsnpfile = ./pindel.out.current_final.gvip.dbsnp_present.vcf\n";
	print PINDELP "EOF\n";
	print PINDELP "myRUNDIR=".$sample_full_path."/pindel\n";
	print PINDELP "pindelout=\${RUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
	print PINDELP "statfile=complete.pindel.parser\n";
	print PINDELP "localstatus=\${myRUNDIR}\/status\/\${statfile}\n";
	print PINDELP "if [ ! -d \${myRUNDIR}\/status ]\n";
	print PINDELP "then\n";
	print PINDELP "mkdir \${myRUNDIR}\/status\n";
	print PINDELP "fi\n";
	print PINDELP "if [ $status_rerun -eq 1 ]\n";
	print PINDELP "  then\n";
	print PINDELP "rm \${localstatus}\n";
	print PINDELP "fi\n";
	print PINDELP "if [ ! -f  \${localstatus} ]\n";
	print PINDELP "then\n";
	print PINDELP "cd \${RUNDIR}/pindel\n";
	print PINDELP "outlist=pindel.out.filelist\n";
	print PINDELP "find \. -name \'*_D\' -o -name \'*_SI\' -o -name \'*_INV\' -o -name \'*_TD\'  > \./\${outlist}\n";
	print PINDELP 'list=$(xargs -a  ./$outlist)'."\n";
	print PINDELP "pin_var_file=pindel.out.raw\n";
	print PINDELP 'cat $list | grep ChrID > ./$pin_var_file'."\n";
	print PINDELP "     ".$run_script_perl."pindel_filter.v0.5.pl ./pindel_filter.input\n"; 
	print PINDELP 'pre_current_final=$pin_var_file.CvgVafStrand_pass.Homopolymer_pass.vcf'."\n";
	print PINDELP 'for mytmp in $pin_var_file.CvgVafStrand_pass.vcf  $pre_current_final  ${pre_current_final/%pass.vcf/fail.vcf} ; do'."\n";
	print PINDELP "     ".$run_script_perl.'genomevip_label.pl Pindel ./$mytmp ./${mytmp/%vcf/gvip.vcf}'."\n";
	print PINDELP "done\n";
	print PINDELP 'current_final=${pin_var_file/%raw/current_final.gvip.Somatic.vcf}'."\n";
	print PINDELP 'cat ./${pre_current_final/%vcf/gvip.vcf} > ./$current_final'."\n";
	print PINDELP "export JAVA_HOME=$java_dir\n";
	print PINDELP "export JAVA_OPTS=\"-Xmx10g\"\n";
	print PINDELP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print PINDELP "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
	print PINDELP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
	print PINDELP "else\n";
	print PINDELP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
	print PINDELP "fi\n";
	print PINDELP "     ".$run_script_perl."dbsnp_filter.pl \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input\n";	
	print PINDELP '		grep "Error occurred during initialization of VM" ${pindelout}',"\n"; 
	print PINDELP '		CHECK=$?',"\n";
	print PINDELP '		while [ ${CHECK} -eq 0 ]',"\n";
	print PINDELP "		do\n";	
	print PINDELP "     ".$run_script_perl."dbsnp_filter.pl \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input\n";
	print PINDELP '      grep "Error occurred during initialization of VM" ${pindelout}',"\n";
	print PINDELP '			CHECK=$?',"\n";
	print PINDELP "		done\n";
	print PINDELP "touch \${localstatus}\n";
	print PINDELP "fi\n";
	close PINDELP;
	
	my $sh_file=$job_files_dir."/".$current_job_file;
	
	### Submit job to background

	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}

##########################################################################################
### Step 7: Run Merge VCF (merge output from three parsing steps and annotate with VEP)
##########################################################################################

sub sub_merge_vcf{

	$current_job_file = "j7_merge_vcf_".$sample_name.".sh";
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

	### Specify error and log files

	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;
	
	my $hg19="Hg19";
	
	### Print out the shell script (job script)
	
	open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
	
	print MERGE "#!/bin/bash\n";
	print MERGE "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
	print MERGE "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
	print MERGE "myRUNDIR=".$sample_full_path."/varscan\n";
	print MERGE "STATUSDIR=".$sample_full_path."/status\n";
	print MERGE "RUNDIR=".$sample_full_path."\n";
	print MERGE "export SAMTOOLS_DIR=".$SAMTOOLS_DIR."\n";
	print MERGE "export JAVA_HOME=$java_dir\n";
	print MERGE "export JAVA_OPTS=\"-Xmx10g\"\n";
	print MERGE "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print MERGE "STRELKA_VCF="."\${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
	print MERGE "VARSCAN_VCF="."\${RUNDIR}/varscan/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print MERGE "PINDEL_VCF="."\${RUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
	print MERGE "VARSCAN_INDEL="."\${RUNDIR}/varscan/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf\n";
	print MERGE "MERGER_OUT="."\${RUNDIR}/merged.vcf\n";
	print MERGE "cat > \${RUNDIR}/vep.merged.input <<EOF\n";
	print MERGE "merged.vep.vcf = ./merged.filtered.vcf\n"; 
	print MERGE "merged.vep.output = ./merged.VEP.vcf\n";
	print MERGE "merged.vep.vep_cmd =".$vepscript."\n";
	print MERGE "merged.vep.cachedir =".$vepcache."\n";
	print MERGE "merged.vep.reffasta = $f_ref_annot\n";
	print MERGE "merged.vep.assembly = GRCh37\n";
	print MERGE "EOF\n";
	print MERGE "java \${JAVA_OPTS} -jar $gatk -R $h37_REF -T CombineVariants -o \${MERGER_OUT} --variant:varscan \${VARSCAN_VCF} --variant:strelka \${STRELKA_VCF} --variant:varindel \${VARSCAN_INDEL} --variant:pindel \${PINDEL_VCF} -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,pindel,varindel\n";
	print MERGE "if [ $ref_name = $hg19 ]\n";
	print MERGE "then\n";	
	print MERGE "     ".$run_script_perl."vaf_filter_v1.1.pl \${RUNDIR}\n";
	print MERGE "else\n";
	print MERGE "     ".$run_script_perl."vaf_filter_v1.1.pl \${RUNDIR}\n";
	print MERGE "fi\n";
	print MERGE "cd \${RUNDIR}\n";
	print MERGE ". $run_script_path/set_envvars\n"; 
	print MERGE "     ".$run_script_perl."vep_annotator.pl ./vep.merged.input >&./vep.merged.log\n";	
	close MERGE;

	my $sh_file=$job_files_dir."/".$current_job_file;

	### Submit job to background
	
	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}

##########################################################################################
### Step 8: Run vcf2maf
##########################################################################################

sub sub_vcf_2_maf{
	
	$current_job_file = "j8_vcf_2_maf_".$sample_name.".sh";
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

	### Specify error and log files

	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;

	### Print out the shell script (job script)

	open(MAF, ">$job_files_dir/$current_job_file") or die $!;
	
	print MAF "#!/bin/bash\n";
	print MAF "F_VCF_1=".$sample_full_path."/merged.filtered.vcf\n";
	print MAF "F_VCF_2=".$sample_full_path."/".$sample_name.".vcf\n";
	print MAF "F_VEP_1=".$sample_full_path."/merged.VEP.vcf\n";
	print MAF "F_VEP_2=".$sample_full_path."/".$sample_name.".vep.vcf\n";
	print MAF "F_maf=".$sample_full_path."/".$sample_name.".maf\n";
	print MAF "rm \${F_VCF_2}\n";
	print MAF "rm \${F_VEP_2}\n"; 
	print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
	print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
	print MAF "     ".$run_script_perl."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf	\${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot --filter-vcf $f_exac\n";
	#print MAF "     ".$run_script_perl."splice_site_check.pl $sample_full_path\n"; 
	close MAF;
	
	my $sh_file=$job_files_dir."/".$current_job_file;

	### Submit job to background

	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}

##########################################################################################
### Step 9: Generate final report
##########################################################################################

sub sub_gen_report{
	
	$current_job_file = "j9_gen_report_".$sample_name.".sh";
	
	### Specify error and log files
	
	my $lsf_out=$log_file_dir."/".$current_job_file.".out";
	my $lsf_err=$log_file_dir."/".$current_job_file.".err";
	`rm $lsf_out`;
	`rm $lsf_err`;
	#`rm $current_job_file`;
	
	### Print out the shell script (job script)
	
	open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
	
	print REPRUN "#!/bin/bash\n";
	print REPRUN "		".$run_script_perl."generate_final_report.pl ".$run_dir."\n";
	close REPRUN;
	
	my $sh_file=$job_files_dir."/".$current_job_file;

	### Submit job to background
	
	$sub_com = "nohup sh $sh_file > $lsf_out 2> $lsf_err &";
	print $sub_com;
	system ($sub_com);
}
