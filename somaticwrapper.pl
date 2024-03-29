######### Song Cao###########
##### email: scao@wustl.edu ####
## pipeline for somatic variant callings ##
#	somatic_variant_callings.pl #
###	updated date: 04/05/2017 ###
### updated date: 04/18/2017 ###
### add vcf2maf.pl ###
### 07/14/2017 ##
##3 vaf_filter.pl ###
### 08/25/2017 ####
### add docker env for mgi ##
### 09/28/17##
## add finishing checks##
### 11/13/17 ###
## add option for log directory##
## 09/26/18 ##
## use mutect1.7; merging using 2 over 3 callers ##
## add tsl for vcf2maf.pl ##

## 08/19/19 ##
## add a checking step for vcf before merging ##

## add smg recovering module, Apr 23, 2020 ##
## cp1 version, Oct 06, 2021

#!/usr/bin/perl
##!/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;

my $version = 2.2;

#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";
#usage information

(my $usage = <<OUT) =~ s/\t+//g;
Somatic variant calling pipeline 
Pipeline version: $version

$yellow     
Usage: perl $0  --srg --sre --wgs --rdir --ref --log --q --mincovt --mincovn --minvaf --maxindsize --exonic --smg --groupname --users --step 

$normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<wgs> = 1 if it is wgs data and otherwise it is 0; If you want to output the maf for all variants, set exonic to 0
<groupname> = job group name
<users> = user name for job group
<srg> = bam having read group or not: 1, yes and 0, no (default 1)
<sre> = re-run and overwrite previous results: 1, yes and 0, no  (default 0)
>>>>>>> a4a79c27f83ba1ebce51a18020e77b941eb33d2b
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
<q> which queue for submitting job; research-hpc, ding-lab, long (default)
<mincovt> minimum coverage for tumor: default >=14
<mincovn> minimum coverage for normal: default >=8
<minvaf> minimum somatic vaf: default >=0.05
<maxindsize> default <=100
<exonic> output exonic region: 1 Yes, 0 No
<smg> use smg list for calling
hg38: /storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
 
$green [1]  Run streka
$green [2]  Run Varscan
$green [3]  Run Pindel
$green [4]  Run mutect
$yellow [5]  Parse mutect result
$yellow [6]  Parse streka result
$yellow [7]  Parse VarScan result
$yellow [8]  Parse Pindel
$cyan [9]  QC vcf files  
$cyan [10] Merge vcf files  
$cyan [11] Generate maf file 
$cyan [12] Generate merged maf file
$cyan [13] Annotate dnp and remove nearby snv near an indel
$red [14] Clean unnecessary intermediate files
$normal

OUT

#__DEFAULT NUMBER OF BINS IE (MUST BE INTEGER)
my $step_number = -1;
my $status_rg = 1;
my $status_rerun=0; 
my $status_exonic=1; 

#__HELP (BOOLEAN, DEFAULTS TO NO-HELP)
my $help = 0;
my $q_name="";
my $s_wgs="";
#__FILE NAME (STRING, NO DEFAULT)
my $run_dir="";
my $log_dir="";
my $h38_REF="";
my $db_smg="";
#my $ref_name="";
my $chr_status=0;
my $maxindsize=100; 
my $mincov_t=14; 
my $mincov_n=8; 
my $minvaf=0.05;
my $compute_username="";
my $group_name="";

#__PARSE COMMAND LINE
my $status = &GetOptions (
      "step=i" => \$step_number,
      "srg=i" => \$status_rg,
      "sre=i" => \$status_rerun,	
      "wgs=i"  => \$s_wgs, 
      "groupname=s" => \$group_name,
      "users=s" => \$compute_username,	
       "exonic=i" => \$status_exonic,
      "rdir=s" => \$run_dir,
	  "ref=s"  => \$h38_REF,
	  "smg=s" => \$db_smg,
	  "log=s"  => \$log_dir,
	  "q=s" => \$q_name,
	  "mincovt=i"  => \$mincov_t,
      "mincovn=i"  => \$mincov_n,		
	  "minvaf=f"  => \$minvaf,
	  "maxindsize=i"  => \$maxindsize,
      "log=s"  => \$log_dir,
      "q=s" => \$q_name,
   	  "help" => \$help, 
	);
 
#print $status,"\n";

if ($help || $run_dir eq "" || $log_dir eq "" || $group_name eq "" || $compute_username eq "" || $step_number<0 || $db_smg eq "") {
	 print "wrong option\n";
	  print $usage;
      exit;
   }

print "run dir=",$run_dir,"\n";
print "log dir=",$log_dir,"\n";
print "step num=",$step_number,"\n";
print "status rerun=",$status_rerun,"\n";
print "status exonic=",$status_exonic,"\n";
print "status readgroup=",$status_rg,"\n";
print "queue name=",$q_name,"\n";
print "minvaf = ",$minvaf,"\n"; 
print "job group=",$group_name,"\n";
print "user group=",$compute_username,"\n";

#<STDIN>; 
#exit;

if($q_name eq "") 
{
	$q_name="long";
}

if($s_wgs eq "") 
{
	$s_wgs=0; 
}

if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

die $usage unless ($step_number >=0)&&(($step_number <= 14));
my $email = "scao\@wustl\.edu";

my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-1];
my $HOME1=$log_dir;

if (! -d $HOME1)
{
`mkdir $HOME1`; 
}
if (! -d $HOME1."/tmpsomatic") {
    `mkdir $HOME1"/tmpsomatic"`;
}
my $job_files_dir = $HOME1."/tmpsomatic";
#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR_SOMATIC") {
    `mkdir $HOME1"/LSF_DIR_SOMATIC"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR_SOMATIC";
#GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip
# obtain script path
#my $script_dir="/gscuser/scao/scripts/git/somaticwrapper";
my $run_script_path =`echo \$PWD`;
chomp $run_script_path;
#my $run_script_path = `dirname $0`;
my $script_dir=$run_script_path; 
print $script_dir,"\n";
#<STDIN>;
#my $run_script_path=$script_dir; 
#chomp $run_script_path;

#$run_script_path = "/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0/bin/perl ".$run_script_path."/";
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";

### running tools: USER needs to change according where the tools are installed.##

#my $mutect="/gscuser/scao/tools/mutect-1.1.7.jar";
my $STRELKA_DIR2="/storage1/fs1/songcao/Active/Software/strelka-2.9.10.centos6_x86_64/bin";
my $pindel="/storage1/fs1/songcao/Active/Software/anaconda3/bin/pindel";
my $PINDEL_DIR="/storage1/fs1/songcao/Active/Software/anaconda3/bin";
my $picardexe="/storage1/fs1/songcao/Active/Software/picard/picard.jar";
my $gatk="/storage1/fs1/songcao/Active/Software/GenomeAnalysis/GenomeAnalysisTK.jar";
my $java_dir="/storage1/fs1/songcao/Active/Software/jre1.8.0_121";
my $java_mutect="/storage1/fs1/songcao/Active/Software/jre1.7.0_80";
my $snpsift="/storage1/fs1/songcao/Active/Software/snpEff/20150522/SnpSift.jar";
my $gatkexe3="/storage1/fs1/songcao/Active/Software/gatk/3.7/GenomeAnalysisTK.jar";
my $mutect1="/storage1/fs1/songcao/Active/Software/mutect/mutect-1.1.7.jar";
my $samtools="/storage1/fs1/songcao/Active/Software/samtools/1.2/bin";
my $samtoolsexe="/storage1/fs1/songcao/Active/Software/samtools/1.2/bin/samtools";
my $varscan="/storage1/fs1/songcao/Active/Software/varscan/2.3.8.ndown";
my $bamreadcount="/storage1/fs1/songcao/Active/Software/bam-readcount/0.7.4/bam-readcount";
#my $vepannot="/storage1/fs1/songcao/Active/Database/hg38_database/vep/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.v102.pl";
my $vepannot="/storage1/fs1/dinglab/Active/Projects/scao/gitshared/ensembl-vep/vep";

#my $vepcache="/storage1/fs1/songcao/Active/Database/hg38_database/vep/v85";

my $vepcache="/storage1/fs1/songcao/Active/Database/hg38_database/vep/v102";

my $DB_SNP_NO_CHR="/storage1/fs1/songcao/Active/Database/hg38_database/DBSNP/00-All.vcf";
my $DB_SNP="/storage1/fs1/songcao/Active/Database/hg38_database/DBSNP/00-All.chr.vcf";

my $DB_COSMIC="/storage1/fs1/songcao/Active/Database/hg38_database/cosmic/CosmicAllMuts.HG38.sort.chr.vcf";
my $DB_COSMIC_NO_CHR="/storage1/fs1/songcao/Active/Database/hg38_database/cosmic/CosmicAllMuts.HG38.sort.vcf"; 

my $DB_SNP_NO_COSMIC="/storage1/fs1/songcao/Active/Database/hg38_database/cosmic/00-All.HG38.pass.cosmic.vcf";
my $f_ref_annot="/storage1/fs1/songcao/Active/Database/hg38_database/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa";
my $TSL_DB="/storage1/fs1/songcao/Active/Database/tsl/wgEncodeGencodeTranscriptionSupportLevelV23.txt";
my $h38_REF_bai=$h38_REF.".fai";
my $f_gtf= "/storage1/fs1/songcao/Active/Database/hg38_database/GTF/Homo_sapiens.GRCh38.85.gtf";

my $first_line=`head -n 1 $h38_REF`;

if($first_line=~/^\>chr/) { $chr_status=1; }

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure
#&check_input_dir($run_dir);
# start data processsing

#`bgadd -L 70 $compute_username/$group_name`;


if (($step_number < 12 && $step_number>0) || $step_number == 14) {
    #begin to process each sample
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
			#	sleep 1;
                $current_job_file="";
                if($step_number==0)
                {  
			&bsub_strelka();
			&bsub_varscan();
			&bsub_pindel();
			&bsub_mutect();
			&bsub_parse_mutect(); 
			&bsub_parse_strelka();
			&bsub_parse_varscan();
			&bsub_parse_pindel();
			&bsub_merge_vcf();
			&bsub_vcf_2_maf();
		}elsif ($step_number == 1) {
                    &bsub_strelka();
                }elsif ($step_number == 2) {
                    &bsub_varscan(1);
                }elsif ($step_number == 3) {
		    &bsub_pindel(1);
                }elsif ($step_number == 4){
                    &bsub_mutect(1);
                }elsif ($step_number == 5){
                    &bsub_parse_mutect(1);
                }elsif ($step_number == 6) {
		    &bsub_parse_strelka(1);
                }elsif ($step_number == 7) {
		    &bsub_parse_varscan(1);
                }elsif ($step_number == 8) {
                    &bsub_parse_pindel(1);
                }elsif ($step_number == 9) {
                    &bsub_qc_vcf(1);
                }elsif ($step_number == 10) {
                    &bsub_merge_vcf(1);
                }elsif ($step_number == 11) {
                    &bsub_vcf_2_maf(1);
                }elsif ($step_number == 14) {
		    &bsub_clean(1);
				} 
           }
        }
    }
}

if($step_number==12)
    {

    print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
    $hold_job_file=$current_job_file; 
    $current_job_file = "j12_Run_report_".$working_name.".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    my $working_name= (split(/\//,$run_dir))[-1];
    my $f_maf=$run_dir."/".$working_name.".withmutect.maf";
    my $f_maf_rc=$f_maf.".rc";
    my $f_maf_rc_caller=$f_maf_rc.".caller";
    open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
    print REPRUN "#!/bin/bash\n";
    print REPRUN "		".$run_script_path."generate_final_report.pl ".$run_dir." ".$status_exonic."\n";
    print REPRUN "      ".$run_script_path."add_rc.pl ".$run_dir." ".$f_maf." ".$f_maf_rc."\n";
    print REPRUN "      ".$run_script_path."add_caller.pl ".$run_dir." ".$f_maf_rc." ".$f_maf_rc_caller."\n";
    close REPRUN;

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);

}

### Annotate dnp 
### keep indel (for cocoexistence of indel and snv) 

print "annotation\n"; 

if($step_number==13)
    {

    print $yellow, "annotate dnp and remove snv near an indel",$normal, "\n";
    $hold_job_file=$current_job_file;
    $current_job_file = "j13_dnp_".$working_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    #`rm $current_job_file`;
    my $working_name= (split(/\//,$run_dir))[-1];
    my $f_maf=$run_dir."/".$working_name.".withmutect.maf.rc.caller";
    my $f_maf_rm_snv=$run_dir."/".$working_name.".remove.nearby.snv.maf";
    my $f_maf_removed=$run_dir."/".$working_name.".remove.nearby.snv.maf.removed";
    my $f_maf_dnp_tmp=$run_dir."/".$working_name.".dnp.annotated.tmp.maf";
    my $f_maf_dnp_tmp_merge=$run_dir."/".$working_name.".dnp.annotated.tmp.maf.merge";
    my $f_maf_dnp=$run_dir."/".$working_name.".dnp.annotated.maf";
    my $f_maf_coding_dnp=$run_dir."/".$working_name.".dnp.annotated.coding.maf";

    my $f_bam_list=$run_dir."/input.bam.list";

    open(OUTB,">$f_bam_list"); 

	foreach my $s (`ls $run_dir`) 
	{

	my $str=$s; 
	chomp($str);
 
	my $dir_2=$run_dir."/".$str; 

	if(-d $dir_2) 
	{
 	my $f_bam=$dir_2."/".$str.".T.bam";

	if(-e $f_bam) { print OUTB $str,"_T","\t",$f_bam,"\n"; }
  				
	}
			
	}
	
    open(DNP, ">$job_files_dir/$current_job_file") or die $!;

    print DNP "#!/bin/bash\n";
    print DNP "      ".$run_script_path."remove_nearby_snv.pl $f_maf $f_maf_rm_snv"."\n";
    print DNP "      ".$run_script_path."cocoon.pl $f_maf_rm_snv $f_maf_dnp_tmp $log_dir --bam $f_bam_list --samt $samtoolsexe --merge --genome $h38_REF --gtf $f_gtf --snvonly"."\n";
    print DNP "		 ".$run_script_path."add_dnp.pl $f_maf_rm_snv $f_maf_dnp_tmp_merge $f_maf_dnp"."\n";

    print DNP "          ".$run_script_path."generate_coding_report.pl $f_maf_dnp $f_maf_coding_dnp"."\n";
### remove tmp files ##
    print DNP "rm $f_maf_dnp_tmp_merge\n";
    print DNP "rm $f_maf_dnp_tmp\n";
    print DNP "rm $f_maf_rm_snv\n"; 
    print DNP "rm $f_maf_removed\n";
	
	#print DNP "      ".$run_script_path."add_caller.pl ".$run_dir." ".$f_maf_rc." ".$f_maf_rc_caller."\n";
    close DNP;

    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ($bsub_com);


    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>100000] rusage[mem=100000]\" -M 100000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;

    system ($bsub_com);

}


exit;

## run strelka

sub bsub_strelka{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_streka_".$sample_name.".sh"; 
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
	my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    if (! -e $IN_bam_T) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_T does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_T) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    }
	if (! -e $IN_bam_N) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_N does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_N) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_N is empty!", $normal, "\n\n";
    }
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
	#`rm $current_job_file`;

    open(STREKA, ">$job_files_dir/$current_job_file") or die $!;
    print STREKA "#!/bin/bash\n";
    #print STREKA "#BSUB -n 1\n";
    #print STREKA "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print STREKA "#BSUB -M 30000000\n";
    #print STREKA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print STREKA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print STREKA "#BSUB -J $current_job_file\n";
    #print STREKA "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print STREKA "#BSUB -q long\n";
    #print STREKA "#BSUB -q research-hpc\n";
	#print STREKA "#BSUB -q long\n";
	#print STREKA "scr_t0=\`date \+\%s\`\n";
    print STREKA "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print STREKA "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print STREKA "myRUNDIR=".$sample_full_path."/strelka\n";
    print STREKA "STATUSDIR=".$sample_full_path."/status\n";
    print STREKA "RESULTSDIR=".$sample_full_path."/results\n"; 
    print STREKA "SG_DIR=".$sample_full_path."/strelka\n"; 
    print STREKA "RUNDIR=".$sample_full_path."\n";
    print STREKA "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
    print STREKA "STRELKA_VCF=".$sample_full_path."/strelka/strelka_out/results/somatic.snvs.vcf.gz"."\n"; 
    #print STREKA "STRELKA_VCF=".$sample_full_path."/strelka/strelka_out/results/passed.somatic.snvs.vcf"."\n";   
	#print STREKA "CONFDIR="."/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/config\n";
    print STREKA "TASK_STATUS=".$sample_full_path."/strelka/strelka_out/task.complete"."\n";
    print STREKA "export SAMTOOLS_DIR=$samtools\n";
    print STREKA "export JAVA_HOME=$java_dir\n";
    print STREKA "export JAVA_OPTS=\"-Xmx10g\"\n";
    print STREKA "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print STREKA "if [ ! -d \${myRUNDIR} ]\n";
    print STREKA "then\n";
    print STREKA "mkdir \${myRUNDIR}\n";
    print STREKA "fi\n";
	### re-run, then delete task.complete file ###
    print STREKA "if [ $status_rerun -eq 1 ]\n";
    print STREKA "  then\n";
    print STREKA "rm \${TASK_STATUS}\n";
    print STREKA "fi\n";
    ## if STRELKA_VCF not existed 
    print STREKA "if [ ! -f \${STRELKA_VCF} ]\n";
    print STREKA "  then\n";
    print STREKA "rm \${TASK_STATUS}\n";
    print STREKA "fi\n";

    print STREKA "if [ ! -f  \${TASK_STATUS} ]\n";
    print STREKA "then\n";
    print STREKA "if [ -d \${STRELKA_OUT} ]\n";
    print STREKA "then\n";
    print STREKA "rm -rf \${STRELKA_OUT}\n";
    print STREKA "fi\n";
    print STREKA "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n"; 
    print STREKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print STREKA "else\n";
    print STREKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print STREKA "fi\n";
	#print STREKA "put_cmd=\"ln -s\"\n";
	#print STREKA "del_cmd=\"rm -f\"\n";
	#print STREKA "del_local=\"rm -f\"\n";
	#print STREKA "statfile=incomplete.strelka\n";
	#print STREKA "localstatus=".$sample_full_path."/status/\$statfile\n";
	#print STREKA "touch \$localstatus\n";
    #print STREKA ". $script_dir/set_envvars\n";
#   	print STREKA ". /gscmnt/gc2525/dinglab/rmashl/Software/perl/set_envvars\n";
   # print STREKA 	
    print STREKA "if [ $s_wgs -eq 1 ]\n";
    print STREKA "  then\n";
    print STREKA "   ".$STRELKA_DIR2."/configureStrelkaSomaticWorkflow.py --normalBam=\$NBAM --tumorBam=\$TBAM --referenceFasta=$h38_REF --callMemMb=1024 --runDir=\$STRELKA_OUT\n";
    print STREKA "else\n";	
    print STREKA "   ".$STRELKA_DIR2."/configureStrelkaSomaticWorkflow.py --normalBam=\$NBAM --tumorBam=\$TBAM --referenceFasta=$h38_REF --callMemMb=1024 --exome --runDir=\$STRELKA_OUT\n";
	#print STREKA "   ".$STRELKA_DIR2."/configureStrelkaWorkflow.pl --normal \$NBAM --tumor \$TBAM --ref ". $h38_REF." --config $script_dir/strelka.ini --output-dir \$STRELKA_OUT\n";
    print STREKA "fi\n";
    print STREKA "cd \$STRELKA_OUT\n";
    print STREKA "\$STRELKA_OUT/runWorkflow.py -m local -j 8 -g 4\n";
    print STREKA "touch \${TASK_STATUS}\n";
    print STREKA "fi\n";
    close STREKA;

    my $sh_file=$job_files_dir."/".$current_job_file;

	if($q_name eq "research-hpc")
	{
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";     }
	else { 
	    $bsub_com = "bsub -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -o $lsf_out -e $lsf_err bash $sh_file\n"; 
	}
    #$bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";     

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

}

## run varscan 

sub bsub_varscan{

	my ($step_by_step) = @_;
	if ($step_by_step) {
		$hold_job_file = "";
	}else{
		$hold_job_file = $current_job_file;
	}

    $current_job_file = "j2_varscan_".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    if (! -e $IN_bam_T) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_T does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_T) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    }
    if (! -e $IN_bam_N) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_N does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_N) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_N is empty!", $normal, "\n\n";
    }
    open(VARSCAN, ">$job_files_dir/$current_job_file") or die $!;
    print VARSCAN "#!/bin/bash\n";
   # print VARSCAN "#BSUB -n 1\n";
   # print VARSCAN "#BSUB -R \"rusage[mem=30000]\"","\n";
   # print VARSCAN "#BSUB -M 30000000\n";
   # print VARSCAN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
   # print VARSCAN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
   # print VARSCAN "#BSUB -J $current_job_file\n";
  #	print VARSCAN "#BSUB -w \"$hold_job_file\"","\n";
   #   print VARSCAN "#BSUB -q long\n";
    #print VARSCAN "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    #print VARSCAN "#BSUB -q research-hpc\n";
    print VARSCAN "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print VARSCAN "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print VARSCAN "myRUNDIR=".$sample_full_path."/varscan\n";
    print VARSCAN "STATUSDIR=".$sample_full_path."/status\n";
    print VARSCAN "RESULTSDIR=".$sample_full_path."/varscan_results\n";
    print VARSCAN "RUNDIR=".$sample_full_path."\n";
    #print VARSCAN "CONFDIR="."/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/config\n";
    #print VARSCAN "GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip\n";
    print VARSCAN "export VARSCAN_DIR=$varscan\n";
    print VARSCAN "export SAMTOOLS_DIR=$samtools\n";
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
    print VARSCAN "varscan_vcf=\${myRUNDIR}\/varscan.out.som_snv.vcf\n";
    print VARSCAN "if [ ! -d \${myRUNDIR}\/status ]\n";
    print VARSCAN "then\n";
    print VARSCAN "mkdir \${myRUNDIR}\/status\n";
    print VARSCAN "fi\n";
  ### re-run, then delete complete.vs_som_snvindel file ###
    print VARSCAN "if [ $status_rerun -eq 1 ]\n";
    print VARSCAN "  then\n";
    print VARSCAN "rm \${localstatus}\n";
    print VARSCAN "fi\n";
   # check if vcf file exists
    print VARSCAN "if [ ! -f  \${varscan_vcf} ]\n";
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
    print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h38_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somatic - \${TMPBASE} --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal $mincov_n --min-coverage-tumor $mincov_t --min-var-freq $minvaf --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp \${snvoutbase} --output-indel \${indeloutbase} &> \${LOG}\n";
    print VARSCAN '      grep "Error occurred during initialization of VM" ${LOG}',"\n";# one possible blast error (see the end of this script). 
    print VARSCAN '      CHECK=$?',"\n";
    print VARSCAN '      while [ ${CHECK} -eq 0 ]',"\n";
    print VARSCAN "      do\n";
    print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h38_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somatic - \${TMPBASE} --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal $mincov_n --min-coverage-tumor $mincov_t --min-var-freq $minvaf --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp \${snvoutbase} --output-indel \${indeloutbase} &> \${LOG}\n";
    print VARSCAN '      grep "Error occurred during initialization of VM" ${LOG}',"\n";
    print VARSCAN '          CHECK=$?',"\n";
    print VARSCAN "      done\n";
    #print VARSCAN "  fi\n";
    print VARSCAN "touch \${localstatus}\n";
    print VARSCAN "fi\n";
    close VARSCAN;	


    my $sh_file=$job_files_dir."/".$current_job_file;

  #  if($q_name eq "research-hpc")
   # {
   # $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";     }
   # else {        $bsub_com = "bsub -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";   }
   # print $bsub_com;
   # system ($bsub_com);
    
     $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

        print $bsub_com;
        system ($bsub_com);

    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ( $bsub_com );

}


sub bsub_parse_strelka{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j6_parse_strelka".$sample_name.".sh";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    open(STREKAP, ">$job_files_dir/$current_job_file") or die $!;

    print STREKAP "#!/bin/bash\n";
    print STREKAP "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print STREKAP "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print STREKAP "myRUNDIR=".$sample_full_path."/strelka\n";
    print STREKAP "STATUSDIR=".$sample_full_path."/status\n";
    print STREKAP "RESULTSDIR=".$sample_full_path."/results\n";
    print STREKAP "SG_DIR=".$sample_full_path."/strelka\n";
    print STREKAP "RUNDIR=".$sample_full_path."\n";
    print STREKAP "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
    print STREKAP "export SAMTOOLS_DIR=$samtools\n";
    print STREKAP "export VARSCAN_DIR=$varscan\n";
    print STREKAP "export JAVA_HOME=$java_dir\n";
    print STREKAP "export JAVA_OPTS=\"-Xmx10g\"\n";
    print STREKAP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print STREKAP "cat > \${myRUNDIR}/strelka_out/results/strelka_dbsnp_filter.snv.input <<EOF\n";
    print STREKAP "streka.dbsnp.snv.annotator = $snpsift\n";
    print STREKAP "streka.dbsnp.snv.db = $DB_SNP_NO_COSMIC\n";
    print STREKAP "streka.dbsnp.snv.rawvcf = ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
    print STREKAP "streka.dbsnp.snv.mode = filter\n";
    print STREKAP "streka.dbsnp.snv.passfile  = ./strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
    print STREKAP "streka.dbsnp.snv.dbsnpfile = ./strelka.somatic.snv.all.gvip.dbsnp_present.vcf\n";
    print STREKAP "EOF\n";
    print STREKAP "cat > \${myRUNDIR}/strelka_out/results/strelka_dbsnp_filter.indel.input <<EOF\n";
    print STREKAP "streka.dbsnp.indel.annotator = $snpsift\n";
    print STREKAP "streka.dbsnp.indel.db = $DB_SNP_NO_COSMIC\n";
    print STREKAP "streka.dbsnp.indel.rawvcf = ./strelka.somatic.indel.strlk_pass.gvip.vcf\n";
    print STREKAP "streka.dbsnp.indel.mode = filter\n";
    print STREKAP "streka.dbsnp.indel.passfile  = ./strelka.somatic.indel.all.gvip.dbsnp_pass.vcf\n";
    print STREKAP "streka.dbsnp.indel.dbsnpfile = ./strelka.somatic.indel.all.gvip.dbsnp_present.vcf\n";
    print STREKAP "EOF\n";
    print STREKAP "FP_BAM=\`awk \'{if(NR==1){print \$1}}\' \${RUNDIR}/varscan/bamfilelist.inp\`\n";
    print STREKAP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_fpfilter.snv.input <<EOF\n";
    print STREKAP "strelka.fpfilter.snv.bam_readcount = $bamreadcount\n";
    print STREKAP "strelka.fpfilter.snv.bam_file = $IN_bam_T\n";
    print STREKAP "strelka.fpfilter.snv.REF = $h38_REF\n"; 
    print STREKAP "strelka.fpfilter.snv.variants_file = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
    print STREKAP "strelka.fpfilter.snv.passfile = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_pass.vcf\n";
    print STREKAP "strelka.fpfilter.snv.failfile = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_fail.vcf\n";
    print STREKAP "strelka.fpfilter.snv.rc_in = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.in.vcf\n";
    print STREKAP "strelka.fpfilter.snv.rc_out = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.out.vcf\n";
    print STREKAP "strelka.fpfilter.snv.fp_out = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp.out.vcf\n";
    print STREKAP "strelka.fpfilter.snv.min_mapping_qual = 0\n";
    print STREKAP "strelka.fpfilter.snv.min_base_qual = 15\n";
    print STREKAP "strelka.fpfilter.snv.min_num_var_supporting_reads = 4\n";
    print STREKAP "strelka.fpfilter.snv.min_var_allele_freq = $minvaf\n";
    print STREKAP "strelka.fpfilter.snv.min_avg_rel_read_position = 0.10\n";
    print STREKAP "strelka.fpfilter.snv.min_avg_rel_dist_to_3prime_end = 0.10\n";
    print STREKAP "strelka.fpfilter.snv.min_var_strandedness = 0.01\n";
    print STREKAP "strelka.fpfilter.snv.min_allele_depth_for_testing_strandedness = 5\n";
    print STREKAP "strelka.fpfilter.snv.min_ref_allele_avg_base_qual = 30\n";
    print STREKAP "strelka.fpfilter.snv.min_var_allele_avg_base_qual = 30\n";
    print STREKAP "strelka.fpfilter.snv.max_rel_read_length_difference = 0.25\n";
    print STREKAP "strelka.fpfilter.snv.max_mismatch_qual_sum_for_var_reads = 150\n";
    print STREKAP "strelka.fpfilter.snv.max_avg_mismatch_qual_sum_difference = 150\n";
    print STREKAP "strelka.fpfilter.snv.min_ref_allele_avg_mapping_qual = 30\n";
    print STREKAP "strelka.fpfilter.snv.min_var_allele_avg_mapping_qual = 30\n";
    print STREKAP "strelka.fpfilter.snv.max_avg_mapping_qual_difference = 50\n";
    print STREKAP "EOF\n";
    print STREKAP "cd \${STRELKA_OUT}/results\n";
    print STREKAP "		".$run_script_path."get_strelka_passed_snv_indel.pl ./variants/somatic.snvs.vcf.gz ./variants/somatic.indels.vcf.gz ./passed.somatic.snvs.vcf ./passed.somatic.indels.vcf\n";
    print STREKAP "     ".$run_script_path."genomevip_label.pl Strelka ./passed.somatic.snvs.vcf ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
    print STREKAP "     ".$run_script_path."genomevip_label.pl Strelka ./passed.somatic.indels.vcf ./strelka.somatic.indel.strlk_pass.gvip.vcf\n";    
    print STREKAP "strekasnvout=\${STRELKA_OUT}/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
    print STREKAP "strekaindelout=\${STRELKA_OUT}/results/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf\n";
    print STREKAP "statfile=complete.streka_parser\n";
    print STREKAP "localstatus=\${STRELKA_OUT}\/\${statfile}\n";
    print STREKAP "if [ $status_rerun -eq 1 ]\n";
    print STREKAP "  then\n";
    print STREKAP "rm \${localstatus}\n";
    print STREKAP "fi\n";
    print STREKAP "if [ ! -f  \${localstatus} ]\n";
    print STREKAP "then\n"; 
    print STREKAP "     ".$run_script_path."dbsnp_filter.pl ./strelka_dbsnp_filter.snv.input\n";
    print STREKAP "     ".$run_script_path."dbsnp_filter.pl ./strelka_dbsnp_filter.indel.input\n";
    print STREKAP "     ".$run_script_path."snv_filter.pl ./strelka_fpfilter.snv.input\n";
    print STREKAP '      grep "Error occurred during initialization of VM" ${strekasnvout}',"\n";
    print STREKAP '      CHECK=$?',"\n";
    print STREKAP '      while [ ${CHECK} -eq 0 ]',"\n";
    print STREKAP "      do\n";
    print STREKAP "     ".$run_script_path."dbsnp_filter.pl ./strelka_dbsnp_filter.snv.input\n";
    print STREKAP "     ".$run_script_path."dbsnp_filter.pl ./strelka_dbsnp_filter.indel.input\n";
    print STREKAP "     ".$run_script_path."snv_filter.pl ./strelka_fpfilter.snv.input\n";  
    print STREKAP '      grep "Error occurred during initialization of VM" ${strekasnvout}',"\n";
    print STREKAP '          CHECK=$?',"\n";
    print STREKAP "      done\n";
    print STREKAP "touch \${localstatus}\n";
    print STREKAP "  fi\n";
    close STREKAP;

    #my $sh_file=$job_files_dir."/".$current_job_file;

    my $sh_file=$job_files_dir."/".$current_job_file;

  #  $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

$bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

        print $bsub_com;
        system ($bsub_com);

}
sub bsub_parse_varscan{

        my ($step_by_step) = @_;
        if ($step_by_step) {
        $hold_job_file = "";
         }else{
        $hold_job_file = $current_job_file;
        }


  	$current_job_file = "j7_parse_varscan".$sample_name.".sh";

        my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
        my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
        `rm $lsf_out`;
        `rm $lsf_err`;

        open(VARSCANP, ">$job_files_dir/$current_job_file") or die $!;

        print VARSCANP "#!/bin/bash\n";
        print VARSCANP "scr_t0=\`date \+\%s\`\n";
	print VARSCANP "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    	print VARSCANP "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
        print VARSCANP "myRUNDIR=".$sample_full_path."/varscan\n";
        print VARSCANP "STATUSDIR=".$sample_full_path."/status\n";
        print VARSCANP "RUNDIR=".$sample_full_path."\n";
        print VARSCANP "export VARSCAN_DIR=$varscan\n"; 
        print VARSCANP "export SAMTOOLS_DIR=$samtools\n";
        print VARSCANP "export JAVA_HOME=$java_dir\n";
        print VARSCANP "export JAVA_OPTS=\"-Xmx10g\"\n";
        print VARSCANP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
        print VARSCANP "cat > \${RUNDIR}/varscan/vs_dbsnp_filter.snv.input <<EOF\n";
	print VARSCANP "varscan.dbsnp.snv.annotator = $snpsift\n";
	print VARSCANP "varscan.dbsnp.snv.db = $DB_SNP_NO_COSMIC\n"; 
	print VARSCANP "varscan.dbsnp.snv.rawvcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.snv.mode = filter\n";
	print VARSCANP "varscan.dbsnp.snv.passfile  = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.snv.dbsnpfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf\n";
	print VARSCANP "EOF\n";
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_dbsnp_filter.indel.input <<EOF\n";
	print VARSCANP "varscan.dbsnp.indel.annotator = $snpsift\n";
	print VARSCANP "varscan.dbsnp.indel.db = $DB_SNP_NO_COSMIC\n";
	print VARSCANP "varscan.dbsnp.indel.rawvcf = ./varscan.out.som_indel.gvip.Somatic.hc.vcf\n";
	print VARSCANP "varscan.dbsnp.indel.mode = filter\n";
	print VARSCANP "varscan.dbsnp.indel.passfile  = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.indel.dbsnpfile = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_present.vcf\n";
	print VARSCANP "EOF\n";
	print VARSCANP "FP_BAM=\`awk \'{if(NR==2){print \$1}}\' \${RUNDIR}/varscan/bamfilelist.inp\`\n";	
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_fpfilter.somatic.snv.input <<EOF\n";
	print VARSCANP "varscan.fpfilter.snv.bam_readcount = $bamreadcount\n";
	print VARSCANP "varscan.fpfilter.snv.bam_file = \${FP_BAM}\n";
	print VARSCANP "varscan.fpfilter.snv.REF = $h38_REF\n";
	print VARSCANP "varscan.fpfilter.snv.variants_file = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.fpfilter.snv.passfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.fp_pass.vcf\n";
	print VARSCANP "varscan.fpfilter.snv.failfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.fp_fail.vcf\n";
	print VARSCANP "varscan.fpfilter.snv.min_mapping_qual = 0\n";
	print VARSCANP "varscan.fpfilter.snv.min_base_qual = 15\n";
	print VARSCANP "varscan.fpfilter.snv.min_num_var_supporting_reads = 4\n";
	print VARSCANP "varscan.fpfilter.snv.min_var_allele_freq = $minvaf\n";
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
        print VARSCANP "localstatus=\${myRUNDIR}\/status\/\${statfile}\n";
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
	print VARSCANP "if [ ! -f \${snvoutbase}.vcf ]","\n";
	print VARSCANP "  then\n";
	print VARSCANP "mv \${snvoutbase} \${snvoutbase}.vcf","\n";
	print VARSCANP "fi\n";
	print VARSCANP "if [ ! -f \${indeloutbase}.vcf ]","\n";
        print VARSCANP " then\n";
	print VARSCANP "mv \${indeloutbase} \${indeloutbase}.vcf","\n";
        print VARSCANP "fi\n";
	print VARSCANP "     ".$run_script_path."genomevip_label.pl VarScan \${snvoutbase}.vcf  \${snvoutbase}.gvip.vcf\n";
	print VARSCANP "     ".$run_script_path."genomevip_label.pl VarScan \${indeloutbase}.vcf \${indeloutbase}.gvip.vcf\n";
	print VARSCANP "echo \'APPLYING PROCESS FILTER TO SOMATIC SNVS:\' &>> \${LOG}\n";
	print VARSCANP "mysnvorig=./\${snvoutbase}.gvip.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar processSomatic \${mysnvorig} --min-tumor-freq $minvaf --max-normal-freq 0.05 --p-value 0.05  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_path."extract_somatic_other.pl <  \${mysnvorig}  > \${mysnvorig/%vcf/other.vcf}\n";
        print VARSCANP "for kk in Somatic Germline LOH ; do\n";
   	print VARSCANP "thisorig=\${mysnvorig/%vcf/\$kk.vcf}\n";
        print VARSCANP "thispass=\${mysnvorig/%vcf/\$kk.hc.vcf}\n";
   	print VARSCANP "thisfail=\${mysnvorig/%vcf/\$kk.lc.vcf}\n";
	print VARSCANP "done\n";
	print VARSCANP "echo \'APPLYING PROCESS FILTER TO SOMATIC INDELS:\' &>> \$LOG\n";
	print VARSCANP "myindelorig=./\$indeloutbase.gvip.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar processSomatic \${myindelorig}   --min-tumor-freq  $minvaf   --max-normal-freq  0.05   --p-value  0.05  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_path."extract_somatic_other.pl <  \${myindelorig}  >  \${myindelorig/%vcf/other.vcf}\n";
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
### remove mincov cut-off ##
        print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somaticFilter  ./\${thissnvorig} --min-reads2  4   --min-strands2  1   --min-avg-qual  20   --min-var-freq  $minvaf --p-value  0.05   --indel-file  ./\${myindelorig} --output-file  ./\${thissnvpass}  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_path."dbsnp_filter.pl  \${RUNDIR}/varscan/vs_dbsnp_filter.snv.input\n";
	print VARSCANP "     ".$run_script_path."dbsnp_filter.pl \${RUNDIR}/varscan/vs_dbsnp_filter.indel.input\n";
        print VARSCANP "touch \${localstatus}\n";
	print VARSCANP "fi\n";
	close VARSCANP; 

 	my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
	print $bsub_com;
        system ($bsub_com);
	}

sub bsub_pindel{
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j3_pindel".$sample_name.".sh";  
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    open(PINDEL, ">$job_files_dir/$current_job_file") or die $!;
    print PINDEL "#!/bin/bash\n";
    print PINDEL "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print PINDEL "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print PINDEL "myRUNDIR=".$sample_full_path."/pindel\n";
    print PINDEL "CONFIG=\${myRUNDIR}"."/".$sample_name.".config\n";
    print PINDEL "statfile=complete.pindel\n";
    print PINDEL "localstatus=\${myRUNDIR}\/status\/\${statfile}\n";
    print PINDEL "pindel_vcf=\${myRUNDIR}\/".$sample_name."_D\n";
    print PINDEL "if [ ! -d \${myRUNDIR} ]\n";
    print PINDEL "then\n";
    print PINDEL "mkdir \${myRUNDIR}\n";
    print PINDEL "fi\n";
    print PINDEL "if [ ! -d \${myRUNDIR}\/status ]\n";
    print PINDEL "then\n";
    print PINDEL "mkdir \${myRUNDIR}\/status\n";
    print PINDEL "fi\n";
    print PINDEL "if [ $status_rerun -eq 1 ]\n";
    print PINDEL "  then\n";
    print PINDEL "rm \${localstatus}\n";
    print PINDEL "fi\n";
    print PINDEL "if [ ! -f  \${pindel_vcf} ]\n";
    print PINDEL "  then\n";
    print PINDEL "rm \${localstatus}\n";
    print PINDEL "fi\n";
    print PINDEL "if [ ! -f  \${localstatus} ]\n";
    print PINDEL "then\n";
    print PINDEL "echo \"$IN_bam_T\t500\t$sample_name.T\" > \${CONFIG}\n";
    print PINDEL "echo \"$IN_bam_N\t500\t$sample_name.N\" >> \${CONFIG}\n";
    print PINDEL "$pindel -T 4 -f $h38_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"." -m 6 -w 1\n";
    print PINDEL "touch \${localstatus}\n"; 
    print PINDEL "fi\n";
    close PINDEL;

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);


	}


sub bsub_parse_pindel {

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j8_parse_pindel".$sample_name.".sh";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    open(PP, ">$job_files_dir/$current_job_file") or die $!;
    print PP "#!/bin/bash\n";
    print PP "RUNDIR=".$sample_full_path."\n";
    print PP "cat > \${RUNDIR}/pindel/pindel_filter.input <<EOF\n";
    print PP "pindel.filter.pindel2vcf = $PINDEL_DIR/pindel2vcf\n";
    print PP "pindel.filter.variants_file = \${RUNDIR}/pindel/pindel.out.raw\n"; 
    print PP "pindel.filter.REF = $h38_REF\n";
    print PP "pindel.filter.date = 000000\n"; 
    print PP "pindel.filter.heterozyg_min_var_allele_freq = 0.2\n";
    print PP "pindel.filter.homozyg_min_var_allele_freq = 0.8\n";
    print PP "pindel.filter.mode = somatic\n";
    print PP "pindel.filter.apply_filter = true\n";
    print PP "pindel.filter.somatic.min_coverages_t = $mincov_t\n";
    print PP "pindel.filter.somatic.min_coverages_n = $mincov_n\n";
    print PP "pindel.filter.somatic.min_var_allele_freq = $minvaf\n";
## changed true to false ##
    print PP "pindel.filter.somatic.require_balanced_reads = \"false\"\n";
    print PP "pindel.filter.somatic.remove_complex_indels = \"true\"\n";
    print PP "pindel.filter.somatic.max_num_homopolymer_repeat_units = 6\n";
    print PP "EOF\n";
    print PP "cat > \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input <<EOF\n";
    print PP "pindel.dbsnp.indel.annotator = $snpsift\n"; 
    print PP "pindel.dbsnp.indel.db = $DB_SNP_NO_COSMIC\n";
    print PP "pindel.dbsnp.indel.rawvcf = ./pindel.out.current_final.gvip.Somatic.vcf\n";
    print PP "pindel.dbsnp.indel.mode = filter\n";
    print PP "pindel.dbsnp.indel.passfile  = ./pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
    print PP "pindel.dbsnp.indel.dbsnpfile = ./pindel.out.current_final.gvip.dbsnp_present.vcf\n";
    print PP "EOF\n";
    print PP "myRUNDIR=".$sample_full_path."/pindel\n";
    print PP "pindelout=\${myRUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.vcf\n"; 
    print PP "statfile=complete.pindel.parser\n";
    print PP "localstatus=\${myRUNDIR}\/status\/\${statfile}\n";
    print PP "if [ ! -d \${myRUNDIR}\/status ]\n";
    print PP "then\n";
    print PP "mkdir \${myRUNDIR}\/status\n";
    print PP "fi\n";
    print PP "if [ $status_rerun -eq 1 ]\n";
    print PP "  then\n";
    print PP "rm \${localstatus}\n";
    print PP "fi\n";
    print PP "if [ ! -f  \${localstatus} ]\n";
    print PP "then\n";
    print PP "cd \${RUNDIR}/pindel\n";
    print PP "outlist=pindel.out.filelist\n";
    print PP "find \. -name \'*_D\' -o -name \'*_SI\' -o -name \'*_INV\' -o -name \'*_TD\'  > \./\${outlist}\n";
    print PP 'list=$(xargs -a  ./$outlist)'."\n";
    print PP "pin_var_file=pindel.out.raw\n";
    print PP 'cat $list | grep ChrID > ./$pin_var_file'."\n";
    print PP "     ".$run_script_path."pindel_filter.v0.5.pl ./pindel_filter.input\n"; 
    print PP 'pre_current_final=$pin_var_file.CvgVafStrand_pass.Homopolymer_pass.vcf'."\n";
    print PP 'for mytmp in $pin_var_file.CvgVafStrand_pass.vcf  $pre_current_final  ${pre_current_final/%pass.vcf/fail.vcf} ; do'."\n";
    print PP "     ".$run_script_path.'genomevip_label.pl Pindel ./$mytmp ./${mytmp/%vcf/gvip.vcf}'."\n";
    print PP "done\n";
    print PP 'current_final=${pin_var_file/%raw/current_final.gvip.Somatic.vcf}'."\n";
    print PP 'cat ./${pre_current_final/%vcf/gvip.vcf} > ./$current_final'."\n";
#	print PP "fi\n";
    print PP "export JAVA_HOME=$java_dir\n";
    print PP "export JAVA_OPTS=\"-Xmx10g\"\n";
    print PP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print PP "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
    print PP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print PP "else\n";
    print PP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print PP "fi\n";
    print PP "     ".$run_script_path."dbsnp_filter.pl \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input\n";	
    print PP '		grep "Error occurred during initialization of VM" ${pindelout}',"\n"; 
    print PP '		CHECK=$?',"\n";
    print PP '		while [ ${CHECK} -eq 0 ]',"\n";
    print PP "		do\n";	
    print PP "     ".$run_script_path."dbsnp_filter.pl \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input\n";
    print PP '      grep "Error occurred during initialization of VM" ${pindelout}',"\n";
    print PP '			CHECK=$?',"\n";
    print PP "		done\n";
    print PP "touch \${localstatus}\n";
    print PP "fi\n";
    close PP;

    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);	
 
 }

## run mutect ##
sub bsub_mutect{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    $current_job_file = "j4_mutect_".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    if (! -e $IN_bam_T) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_T does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }

    if (! -s $IN_bam_T) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    }

    if (! -e $IN_bam_N) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_N does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }

    if (! -s $IN_bam_N) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_N is empty!", $normal, "\n\n";
    }

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
	`rm $lsf_err`; 

    open(MUTECT, ">$job_files_dir/$current_job_file") or die $!;
    print MUTECT "#!/bin/bash\n";
    print MUTECT "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print MUTECT "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print MUTECT "TBAM_rg=".$sample_full_path."/".$sample_name.".T.rg.bam\n";
    print MUTECT "NBAM_rg=".$sample_full_path."/".$sample_name.".N.rg.bam\n";
    print MUTECT "TBAM_rg_bai=".$sample_full_path."/".$sample_name.".T.rg.bam.bai\n";
    print MUTECT "NBAM_rg_bai=".$sample_full_path."/".$sample_name.".N.rg.bam.bai\n";
    print MUTECT "fcov=".$sample_full_path."/mutect1/mutect.coverage\n";
    print MUTECT "fstat=".$sample_full_path."/mutect1/mutect.status\n";
    print MUTECT "myRUNDIR=".$sample_full_path."/mutect1\n";
    print MUTECT "filtervcf=".$sample_full_path."/mutect1/mutect.filtered.vcf\n";
    print MUTECT "RUNDIR=".$sample_full_path."\n";
    print MUTECT "export SAMTOOLS_DIR=$samtools\n";
    print MUTECT "export JAVA_HOME=$java_mutect\n";
    print MUTECT "export JAVA_OPTS=\"-Xmx10g\"\n";
    print MUTECT "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";

### re-run delete mutect folder and re-run it ##

    print MUTECT "if [ $status_rerun -eq 1 ]\n";
    print MUTECT "  then\n";
    print MUTECT "rm -rf \${myRUNDIR}\n";
    print MUTECT "fi\n";

    print MUTECT "if [ ! -d \${myRUNDIR} ]\n";
    print MUTECT "then\n";
    print MUTECT "mkdir \${myRUNDIR}\n";
    print MUTECT "fi\n";
    print MUTECT "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print MUTECT "if [ $status_rg -eq 0 ]\n";
    print MUTECT "then\n";	
    print MUTECT "java  \${JAVA_OPTS} -jar "."$picardexe AddOrReplaceReadGroups I=\${NBAM} O=\${NBAM_rg} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample_name\n";
    print MUTECT "samtools index \${NBAM_rg}\n";
    print MUTECT "java  \${JAVA_OPTS} -jar "."$picardexe AddOrReplaceReadGroups I=\${TBAM} O=\${TBAM_rg} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample_name\n";
    print MUTECT "samtools index \${TBAM_rg}\n";

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
	my $DB_SNP_MUTECT=$DB_SNP_NO_CHR;
	my $DB_COSMIC_MUTECT=$DB_COSMIC_NO_CHR; 
    	if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
		$DB_SNP_MUTECT=$DB_SNP;
		$DB_COSMIC_MUTECT=$DB_COSMIC; 
		}
        print MUTECT "rawvcf=".$sample_full_path."/mutect1/mutect.raw.$chr.vcf\n";
    	print MUTECT "if [ $status_rerun -eq 1 ]\n";
    	print MUTECT "  then\n";
    	print MUTECT "rm \${rawvcf}\n";
    	print MUTECT "fi\n";
        print MUTECT '  if [ ! -s $rawvcf ]',"\n";
    	print MUTECT "  then\n";
    	print MUTECT "java  \${JAVA_OPTS} -jar "."$mutect1  -T MuTect -R $h38_REF -L $chr1 --dbsnp $DB_SNP_MUTECT --cosmic $DB_COSMIC_MUTECT -I:normal \${NBAM_rg} -I:tumor \${TBAM_rg} --artifact_detection_mode -vcf \${rawvcf}\n";
	    print MUTECT "  fi\n";
	} 

    print MUTECT "rm \${NBAM_rg}\n";
    print MUTECT "rm \${NBAM_rg_bai}\n";
    print MUTECT "rm \${TBAM_rg}\n";
    print MUTECT "rm \${TBAM_rg_bai}\n";
    print MUTECT "else\n";
    foreach my $chr (@chrlist)
    {
	 my $chr1=$chr;
	 my $DB_SNP_MUTECT=$DB_SNP_NO_CHR;
     	 my $DB_COSMIC_MUTECT=$DB_COSMIC_NO_CHR;
         if($chr_status==1) 
         { 
         $chr1="chr".$chr;
         $DB_SNP_MUTECT=$DB_SNP;
         $DB_COSMIC_MUTECT=$DB_COSMIC;
         }


	 print MUTECT "rawvcf=".$sample_full_path."/mutect1/mutect.raw.$chr.vcf\n";	
         print MUTECT "if [ $status_rerun -eq 1 ]\n";
         print MUTECT "  then\n";
         print MUTECT "rm \${rawvcf}\n";
         print MUTECT "fi\n";
	 print MUTECT '  if [ ! -s $rawvcf ]',"\n"; 
         print MUTECT "  then\n";
         print MUTECT "java  \${JAVA_OPTS} -jar "."$mutect1  -T MuTect -R $h38_REF -L $chr1 --dbsnp $DB_SNP_MUTECT --cosmic $DB_COSMIC_MUTECT -I:normal \${NBAM} -I:tumor \${TBAM} --artifact_detection_mode -vcf \${rawvcf}\n";
 	 print MUTECT "  fi\n"; 
      }

      print MUTECT "fi\n";

  
      close MUTECT;

      my $sh_file=$job_files_dir."/".$current_job_file;

      $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
     #my $sh_file=$job_files_dir."/".$current_job_file;
    #if($q_name eq "research-hpc")
    #{
   # $bsub_com = "bsub -q research-hpc -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";     }    
#	else 
#	{        
#	$bsub_com = "bsub -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";                                
 #   }
      print $bsub_com;
      system ($bsub_com);

}


sub bsub_parse_mutect{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j5_parse_mutect.".$sample_name.".sh";

    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    open(PM, ">$job_files_dir/$current_job_file") or die $!;
    print PM "#!/bin/bash\n";
    print PM "export JAVA_HOME=$java_mutect\n";
    print PM "export JAVA_OPTS=\"-Xmx10g\"\n";
    print PM "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print PM "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";  

    foreach my $chr (@chrlist)
    {
    print PM "rawvcf=".$sample_full_path."/mutect1/mutect.raw.$chr.vcf\n";
    print PM "filtervcf=".$sample_full_path."/mutect1/mutect.raw.filtered.$chr.vcf\n";
    print PM "filtervcfsnv=".$sample_full_path."/mutect1/mutect.filter.snv.$chr.vcf\n";
    print PM "filtervcfindel=".$sample_full_path."/mutect1/mutect.filter.indel.$chr.vcf\n";
    print PM "     ".$run_script_path."filter_mutect1.8.pl $samtools/samtools \${rawvcf} \${filtervcf} $mincov_t $mincov_n $minvaf\n";
    print PM "java \${JAVA_OPTS} -jar "."$mutect1  -T SelectVariants -R $h38_REF -V \${filtervcf}  -o  \${filtervcfsnv}  -selectType SNP -selectType MNP"."\n";
    print PM "java \${JAVA_OPTS} -jar "."$mutect1  -T SelectVariants -R $h38_REF -V  \${filtervcf}  -o  \${filtervcfindel}  -selectType INDEL"."\n";
    }

    print PM "     ".$run_script_path."merge_mutect.pl $sample_full_path\n";   
    close PM;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);


}

sub bsub_qc_vcf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j9_qc_vcf.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
	### QC VCF file ##
    open(QC, ">$job_files_dir/$current_job_file") or die $!;
    print QC "#!/bin/bash\n";
    print QC "RUNDIR=".$sample_full_path."\n";	
    print QC "PINDEL_VCF="."\${RUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
    print QC "PINDEL_VCF_QC="."\${RUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.qced.vcf\n";
    print QC "     ".$run_script_path."qc_vcf.pl \${PINDEL_VCF} \${PINDEL_VCF_QC}\n";
    close QC;
    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>100000] rusage[mem=100000]\" -M 100000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);

}


sub bsub_merge_vcf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j10_merge_vcf.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    my $hg19="Hg19"; 
    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    print MERGE "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print MERGE "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print MERGE "myRUNDIR=".$sample_full_path."/varscan\n";
    print MERGE "STATUSDIR=".$sample_full_path."/status\n";
    print MERGE "RUNDIR=".$sample_full_path."\n";
    print MERGE "export SAMTOOLS_DIR=$samtools\n";
    print MERGE "export JAVA_HOME=$java_dir\n";
    print MERGE "export JAVA_OPTS=\"-Xmx10g\"\n";
    print MERGE "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print MERGE "STRELKA_VCF="."\${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
    print MERGE "VARSCAN_VCF="."\${RUNDIR}/varscan/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
    print MERGE "PINDEL_VCF="."\${RUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.qced.vcf\n";
    print MERGE "VARSCAN_INDEL="."\${RUNDIR}/varscan/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf\n";
    print MERGE	"MUTECT_VCF="."\${RUNDIR}/mutect1/mutect.filter.snv.vcf\n";
    print MERGE "STRELKA_INDEL="."\${RUNDIR}/strelka/strelka_out/results/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf\n";; 
    print MERGE "MERGER_OUT="."\${RUNDIR}/merged.withmutect.vcf\n";
    print MERGE "PINDEL_VCF_FILTER="."\${RUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.filtered.vcf\n";
    print MERGE "     ".$run_script_path."filter_large_indel.pl \${PINDEL_VCF} \${PINDEL_VCF_FILTER} $maxindsize\n";
    print MERGE "java \${JAVA_OPTS} -jar $gatk -R $h38_REF -T CombineVariants -o \${MERGER_OUT} --variant:varscan \${VARSCAN_VCF} --variant:strelka \${STRELKA_VCF} --variant:mutect \${MUTECT_VCF} --variant:varindel \${VARSCAN_INDEL} --variant:sindel \${STRELKA_INDEL} --variant:pindel \${PINDEL_VCF_FILTER} -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,mutect,sindel,varindel,pindel\n"; 
    close MERGE;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>100000] rusage[mem=100000]\" -M 100000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);


   }

sub bsub_vcf_2_maf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j11_vcf_2_maf.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    open(MAF, ">$job_files_dir/$current_job_file") or die $!;
    print MAF "#!/bin/bash\n";
  	
    print MAF "F_VCF_1=".$sample_full_path."/merged.withmutect.vcf\n";
    print MAF "F_VCF_rm=".$sample_full_path."/merged.withmutect.rm.largeindel.vcf\n";
    print MAF "F_VCF_1_filtered=".$sample_full_path."/merged.filtered.withmutect.vcf\n";
    print MAF "F_VCF_2=".$sample_full_path."/".$sample_name.".withmutect.vcf\n";
    print MAF "F_VCF_2_filtered=".$sample_full_path."/".$sample_name.".withmutect.filtered.vcf\n";
    print MAF "F_VEP_1=".$sample_full_path."/merged.VEP.withmutect.vcf\n";
    print MAF "F_VEP_1_filtered=".$sample_full_path."/merged.VEP.withmutect.filtered.vcf\n";
    print MAF "F_VEP_2=".$sample_full_path."/".$sample_name.".withmutect.vep.vcf\n";
    print MAF "F_VEP_2_filtered=".$sample_full_path."/".$sample_name.".withmutect.filtered.vep.vcf\n";
    print MAF "F_maf=".$sample_full_path."/".$sample_name.".withmutect.maf\n";
    print MAF "F_maf_filtered=".$sample_full_path."/".$sample_name.".withmutect.filtered.maf\n";
    print MAF "RUNDIR=".$sample_full_path."\n";
	
    print MAF "F_log=".$sample_full_path."/vep.merged.withmutect.log"."\n";
    print MAF "cat > \${RUNDIR}/vep.merged.withmutect.input <<EOF\n";
    print MAF "merged.vep.vcf = ./merged.withmutect.rm.largeindel.vcf\n";
    print MAF "merged.vep.output = ./merged.VEP.withmutect.vcf\n";
    print MAF "merged.vep.vep_cmd = $vepannot\n";
    print MAF "merged.vep.cachedir = $vepcache\n";
    print MAF "merged.vep.reffasta = $f_ref_annot\n";
    print MAF "merged.vep.assembly = GRCh38\n";
    print MAF "EOF\n";
    print MAF "rm \${F_log}\n";
	
    print MAF "F_log_filtered=".$sample_full_path."/vep.merged.withmutect.filtered.log"."\n";
    print MAF "cat > \${RUNDIR}/vep.merged.withmutect.filtered.input <<EOF\n";
    print MAF "merged.vep.vcf = ./merged.filtered.withmutect.vcf\n";
    print MAF "merged.vep.output = ./merged.VEP.withmutect.filtered.vcf\n";
    print MAF "merged.vep.vep_cmd = $vepannot\n";
    print MAF "merged.vep.cachedir = $vepcache\n";
    print MAF "merged.vep.reffasta = $f_ref_annot\n";
    print MAF "merged.vep.assembly = GRCh38\n";
    print MAF "EOF\n";
    print MAF "rm \${F_log_filtered}\n";	
  
	### vep and vcf2maf annotation for all variants to get the annotated gene name for each variant ##
    print MAF "cd \${RUNDIR}\n";
    #print MAF ". $script_dir/set_envvars\n";
    print MAF "     ".$run_script_path."remove_largeindel.pl \${F_VCF_1} \${F_VCF_rm}\n";
    print MAF "     ".$run_script_path."vep_annotator.pl ./vep.merged.withmutect.input >&./vep.merged.withmutect.log\n";
    print MAF "rm \${F_VCF_2}\n";
    print MAF "rm \${F_VEP_2}\n";
    print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
    print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
    print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot --file-tsl $TSL_DB\n";	
	## do the filtering for variants and ignore tumor vaf > 0.05 for gene in smg ##
    print MAF "     ".$run_script_path."vaf_filter_v1.4.pl \${RUNDIR} $sample_name $minvaf $mincov_t $mincov_n $maxindsize $db_smg\n"; 
    print MAF "     ".$run_script_path."vep_annotator.pl ./vep.merged.withmutect.filtered.input >&./vep.merged.withmutect.filtered.log\n";
    print MAF "rm \${F_VCF_2_filtered}\n";
    print MAF "rm \${F_VEP_2_filtered}\n";
    print MAF "ln -s \${F_VCF_1_filtered} \${F_VCF_2_filtered}\n";
    print MAF "ln -s \${F_VEP_1_filtered} \${F_VEP_2_filtered}\n";
    print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2_filtered} --output-maf \${F_maf_filtered} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot --file-tsl $TSL_DB\n"; 
    close MAF;


    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(ensemblorg/ensembl-vep:release_102.0)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);


 }


sub bsub_clean{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j14_clean_".$sample_name.".sh";
    
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    #`rm $current_job_file`;

    open(CLEAN, ">$job_files_dir/$current_job_file") or die $!;
    print CLEAN "#!/bin/bash\n";
    print CLEAN "RP=".$sample_full_path."/pindel/".$sample_name."_RP\n";
    print CLEAN "D=".$sample_full_path."/pindel/".$sample_name."_D\n";
    print CLEAN "TD=".$sample_full_path."/pindel/".$sample_name."_TD\n";
    print CLEAN "INV=".$sample_full_path."/pindel/".$sample_name."_INV\n";
    print CLEAN "SI=".$sample_full_path."/pindel/".$sample_name."_SI\n";	
    print CLEAN "RAW=".$sample_full_path."/pindel/pindel.out.raw\n";  	
    print CLEAN "FAIL=".$sample_full_path."/pindel/pindel.out.raw.CvgVafStrand_fail\n";
    print CLEAN "if [ -f \${RP} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${RP}\n";
    print CLEAN "fi\n";

    print CLEAN "if [ -f \${D} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${D}\n";
    print CLEAN "fi\n";
    
    print CLEAN "if [ -f \${TD} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${TD}\n";
    print CLEAN "fi\n";
    

    print CLEAN "if [ -f \${INV} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${INV}\n";
    print CLEAN "fi\n";
    

    print CLEAN "if [ -f \${SI} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${SI}\n";
    print CLEAN "fi\n";
    
    print CLEAN "if [ -f \${RAW} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${RAW}\n";
    print CLEAN "fi\n";

    print CLEAN "if [ -f \${FAIL} ]\n";
    print CLEAN "then\n";
    print CLEAN "gzip \${FAIL}\n";
    print CLEAN "fi\n";

    close CLEAN;


    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>100000] rusage[mem=100000]\" -M 100000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);

}
 
