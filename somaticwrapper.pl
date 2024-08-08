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

my $version = 1.6.1;

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
 
$cyan [1] Generate maf file 
$cyan [2] Generate merged maf file
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
my $varscan="/storage1/fs1/songcao/Active/Software/varscan/2.3.8.ndown";
my $bamreadcount="/storage1/fs1/songcao/Active/Software/bam-readcount/0.7.4/bam-readcount";
#my $vepannot="/storage1/fs1/songcao/Active/Database/hg38_database/vep/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl";
#my $vepcache="/storage1/fs1/songcao/Active/Database/hg38_database/vep/v85";
my $vepannot="/storage1/fs1/dinglab/Active/Projects/scao/gitshared/ensembl-vep/vep";
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


if ($step_number == 1) {
    #begin to process each sample
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
			#	sleep 1;
                $current_job_file="";
                if($step_number==1)
                {  
	          &bsub_vcf_2_maf();
	        } 
           }
        }
    }
}

if($step_number == 2)
    {

	print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
	$hold_job_file=$current_job_file; 
	$current_job_file = "j2_Run_report_".$working_name.".sh"; 
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
#    print REPRUN "      ".$run_script_path."add_rc.pl ".$run_dir." ".$f_maf." ".$f_maf_rc."\n";
 #   print REPRUN "      ".$run_script_path."add_caller.pl ".$run_dir." ".$f_maf_rc." ".$f_maf_rc_caller."\n";
    close REPRUN;

     my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);

}



exit;


sub bsub_vcf_2_maf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j1_vcf_2_maf.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    open(MAF, ">$job_files_dir/$current_job_file") or die $!;
    print MAF "#!/bin/bash\n";
    print MAF "F_VCF_1_gz=".$sample_full_path."/merged.withmutect.vcf.gz\n";
    print MAF "F_VCF_1=".$sample_full_path."/merged.withmutect.vcf\n";
    print MAF "F_vcf=".$sample_full_path."/".$sample_name.".vcf\n";
    print MAF "F_vcf_vep=".$sample_full_path."/".$sample_name.".vep.vcf\n";
    print MAF "F_VCF_1_vep=".$sample_full_path."/merged.VEP.withmutect.vcf\n";
    print MAF "F_maf=".$sample_full_path."/".$sample_name.".maf\n";
    print MAF "F_log=".$sample_full_path."/vep.merged.withmutect.log"."\n";
    print MAF "RUNDIR=".$sample_full_path."\n";
    print MAF "cat > \${RUNDIR}/vep.merged.withmutect.input <<EOF\n";
    print MAF "merged.vep.vcf = ./merged.withmutect.vcf\n";
    print MAF "merged.vep.output = ./merged.VEP.withmutect.vcf\n";
    print MAF "merged.vep.vep_cmd = $vepannot\n";
    print MAF "merged.vep.cachedir = $vepcache\n";
    print MAF "merged.vep.reffasta = $f_ref_annot\n";
    print MAF "merged.vep.assembly = GRCh38\n";
    print MAF "EOF\n";
    print MAF "rm \${F_log}\n";
    print MAF "if [ -e \${F_VCF_1_gz} ]\n";
    print MAF "then\n";
    print MAF "gunzip -c \${F_VCF_1_gz} > \${F_VCF_1}\n";
    print MAF "fi\n";
    print MAF "cd \${RUNDIR}\n";
    print MAF "     ".$run_script_path."vep_annotator_all.pl ./vep.merged.withmutect.input >&./vep.merged.withmutect.log\n";
    print MAF "rm \${F_vcf}\n";
    print MAF "rm \${F_vcf_vep}\n";
    print MAF "ln -s \${F_VCF_1} \${F_vcf}\n";
    print MAF "ln -s \${F_VCF_1_vep} \${F_vcf_vep}\n";
    print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_vcf} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot --file-tsl $TSL_DB\n"; 
    close MAF;

    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(ensemblorg/ensembl-vep:release_102.0)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

 }

