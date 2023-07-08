######### Song Cao###########
##### email: scao@wustl.edu ####
## pipeline for tumor-only somatic variant callings ##

#!/usr/bin/perl
##!/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;

my $version = 1.0;

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

$yellow     Usage: perl $0 --chr --step --sre --rdir --ref --log --q --exonic --groupname --users --step 

$normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<srg> = bam having read group or not: 1, yes and 0, no (default 1)
<chr> = with chr or not chr in the reference file: 1, yes and 0, no (default 1) 
<sre> = re-run: 1, yes and 0, no  (default 0)
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
<q> which queue for submitting job; general (default), dinglab
<groupname> which job group to use when submitting jobs on compute1
<users> username on compute1 to use for submitting jobs to job groups
<exonic> output exonic region: 1 Yes, 0 No
 
$green [0] Generate bams if input files are fastqs
$green [1] Run Mutect2
$green [2] Run filter Mutect2 result
$green [3] Run parse Mutect2 result
$cyan [4] Generate maf file 
$cyan [5] Generate merged maf file
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
my $chr_status=1;
my $maxindsize=100; 
my $mincov_t=14; 
#my $mincov_n=8; 
my $minvaf=0.05;
my $compute_username="";
my $group_name="";
#my $compute_group;

#__PARSE COMMAND LINE
my $status = &GetOptions (
      "step=i" => \$step_number,
      "chr=i" => \$chr_status,
      "sre=i" => \$status_rerun,
      "groupname=s" => \$group_name,
      "users=s" => \$compute_username,		
	  "exonic=i" => \$status_exonic,
 #     "group=s" => \$compute_group,
      "rdir=s" => \$run_dir,
	  "ref=s"  => \$h38_REF,
#	  "smg=s" => \$db_smg,
	  "log=s"  => \$log_dir,
	  "q=s" => \$q_name,
#	  "mincovt=i"  => \$mincov_t,
#	  "minvaf=f"  => \$minvaf,
#	  "maxindsize=i"  => \$maxindsize,
      "log=s"  => \$log_dir,
      "q=s" => \$q_name,
   	  "help" => \$help, 
	);
 
#print $status,"\n";

if ($help || $run_dir eq "" || $log_dir eq "" || $step_number<0 || $group_name eq "" || $compute_username eq "") {
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
#print "minvaf = ",$minvaf,"\n"; 

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

#<STDIN>;
#die $usage unless @ARGV == 3;
#my ($run_dir, $status_rg, $step_number) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

die $usage unless ($step_number >=0)&&(($step_number <= 5));
my $email = "scao\@wustl\.edu";
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-1];
my $HOME1=$log_dir;
#store job files here
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
my $run_py_script_path = "python ".$run_script_path."/";

print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";

### running tools: USER needs to change according where the tools are installed.##

my $MINLEN=50;
my $MIN_VAF=0.02;
my $MIN_DEP=8;
my $MIN_MUT=3;

#my $FASTQC="/gscmnt/gc3021/dinglab/hsun/software/miniconda2/bin/fastqc";
#my $Z7="/gscmnt/gc3021/dinglab/hsun/software/miniconda2/bin/7z"; 
# set trimGalore
#my $TRIMGALORE="/gscmnt/gc3021/dinglab/hsun/software/miniconda2/bin/trim_galore";

#TRIMGALORE=/gscmnt/gc3021/dinglab/hsun/software/miniconda2/bin/trim_galore


##================================= Human Genomics
# set for target coverage
#REF=/gscmnt/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa
my $f_known_site="/storage1/fs1/dinglab/Active/Projects/estorrs/pecgs_resources/dnaseq_alignment/dbsnp/00-All.chr.vcf.gz"; 

my $TARGET="/storage1/fs1/songcao/Active/Database/tonlydb/proteinCoding.cds.merged.gencode_hg38_v29.interval_list";
my $GNOMAD_VCF="/storage1/fs1/songcao/Active/Database/tonlydb/af-only-gnomad.hg38.vcf.gz";
# pon normal
my $PANEL_OF_NORMALS_VCF="/storage1/fs1/dinglab/Active/Projects/austins2/tools/somaticwrapper_tonly.v1.0/db/gatk/GDC-gatk4_panel_of_normals/6c4c4a48-3589-4fc0-b1fd-ce56e88c06e4/gatk4_mutect2_4136_pon.hg38.vcf.gz";
# common_biallelic
my $COMMON_BIALLELIC="/storage1/fs1/dinglab/Active/Projects/austins2/tools/somaticwrapper_tonly.v1.0/db/gatk/mutect2_gatk-best-practices.broadIns/af-only-gnomad.hg38.common_biallelic.chr1-22XY.vcf";
# gencode
my $GENOME="/storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa";
#my $GENOME="/storage1/fs1/dinglab/Active/Projects/austins2/tools/somaticwrapper_tonly.v1.0/db/gencode_GRCh38_v29/GRCh38.primary_assembly.genome.fa";
#GENOME=/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data/image.data/A_Reference/GRCh38.d1.vd1.fa
my $GATK="/storage1/fs1/dinglab/Active/Projects/austins2/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar";
#my $JAVA=/gscmnt/gc2737/ding/hsun/software/jre1.8.0_152/bin/java
my $PICARD="/storage1/fs1/songcao/Active/Software/picard/picard.jar";
#not installed
#my $BWA="/gscmnt/gc2737/ding/hsun/software/bwa-0.7.17/bwa";

my $BWA="/storage1/fs1/songcao/Active/Software/anaconda3/bin/bwa"; 
my $SAMTOOLS="/storage1/fs1/songcao/Active/Software/anaconda3/bin/samtools";
my $mutect="/storage1/fs1/songcao/Active/Software/mutect/mutect-1.1.7.jar";
my $STRELKA_DIR2="/storage1/fs1/songcao/Active/Software/strelka-2.9.10.centos6_x86_64/bin";
my $pindel="/storage1/fs1/songcao/Active/Software/anaconda3/bin/pindel";
my $PINDEL_DIR="/storage1/fs1/songcao/Active/Software/anaconda3/bin";
my $picardexe="/storage1/fs1/songcao/Active/Software/picard/picard.jar";
my $gatk="/storage1/fs1/songcao/Active/Software/GenomeAnalysis/GenomeAnalysisTK.jar";
my $java_dir="/storage1/fs1/songcao/Active/Software/jre1.8.0_121";
my $java_bin="/storage1/fs1/songcao/Active/Software/jre1.8.0_121/bin/java";
my $java_mutect="/storage1/fs1/songcao/Active/Software/jre1.7.0_80";
my $snpsift="/storage1/fs1/songcao/Active/Software/snpEff/20150522/SnpSift.jar";
my $gatkexe3="/storage1/fs1/songcao/Active/Software/gatk/3.7/GenomeAnalysisTK.jar";
my $mutect1="/storage1/fs1/songcao/Active/Software/mutect/mutect-1.1.7.jar";
my $samtools="/storage1/fs1/songcao/Active/Software/samtools/1.2/bin";
my $samtoolexe="/storage1/fs1/songcao/Active/Software/samtools/1.2/bin/samtools";
my $varscan="/storage1/fs1/songcao/Active/Software/varscan/2.3.8.ndown";
my $bamreadcount="/storage1/fs1/songcao/Active/bam-readcount/0.7.4/bam-readcount";
my $vepannot="/storage1/fs1/songcao/Active/Database/hg38_database/vep/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl";
#my $vepannot="/storage1/fs1/dinglab/Active/Projects/scao/gitshared/ensembl-vep/vep";
my $vepcache="/storage1/fs1/songcao/Active/Database/hg38_database/vep/v85";
#my $vepcache="/storage1/fs1/songcao/Active/Database/hg38_database/vep/v102";

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

if ($step_number < 5) {
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
		    &bsub_fq2bam();
		} elsif ($step_number == 1) {
                    &bsub_mutect2(1);
                } elsif ($step_number == 2) {
					&bsub_filter_mutect2(1);
                } elsif ($step_number == 3){
                    &bsub_parse_mutect2(1);
                } elsif ($step_number == 4) {
                    &bsub_vcf_2_maf(1);
                } 
           }
        }
    }
}

if($step_number==5)
    {

    print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
    $hold_job_file=$current_job_file; 
    $current_job_file = "j5_Run_report_".$working_name.".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    my $working_name= (split(/\//,$run_dir))[-1];
    my $f_maf=$run_dir."/".$working_name.".mutect2.maf";
    my $f_maf_rc=$f_maf.".rc";
    open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
    print REPRUN "#!/bin/bash\n";
    #print REPRUN "#BSUB -n 1\n";
    #print REPRUN "#BSUB -R \"rusage[mem=40000]\"","\n";
    #print REPRUN "#BSUB -M 40000000\n";
    #print REPRUN "#BSUB -a \'docker(scao/dailybox)\'\n";
    #print REPRUN "#BSUB -q research-hpc\n";
   # print REPRUN "#BSUB -q ding-lab\n";
	#print REPRUN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print REPRUN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print REPRUN "#BSUB -J $current_job_file\n";
	#print REPRUN "#BSUB -w \"$hold_job_file\"","\n";
    print REPRUN "		".$run_script_path."generate_final_report.pl ".$run_dir." ".$status_exonic."\n";
    print REPRUN "      ".$run_script_path."add_rc.pl ".$run_dir." ".$f_maf." ".$f_maf_rc."\n";
    #print REPRUN "      ".$run_script_path."add_caller.pl ".$run_dir." ".$f_maf_rc." ".$f_maf_rc_caller."\n";
    close REPRUN;
    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
	#system ($bsub_com);
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);
}

### Annotate dnp 
### keep indel (for cocoexistence of indel and snv) 

print "annotation\n"; 

if($step_number==6)
    {

    print $yellow, "annotate dnp and remove snv near an indel",$normal, "\n";
    $hold_job_file=$current_job_file;
    $current_job_file = "j6_dnp_".$working_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    #`rm $current_job_file`;
    my $working_name= (split(/\//,$run_dir))[-1];
    my $f_maf=$run_dir."/".$working_name.".mutect2.maf.rc";
    my $f_maf_rm_snv=$run_dir."/".$working_name.".remove.nearby.snv.maf";
    my $f_maf_removed=$run_dir."/".$working_name.".remove.nearby.snv.maf.removed";
    my $f_maf_dnp_tmp=$run_dir."/".$working_name.".dnp.annotated.tmp.maf";
    my $f_maf_dnp_tmp_merge=$run_dir."/".$working_name.".dnp.annotated.tmp.maf.merge";
    my $f_maf_dnp=$run_dir."/".$working_name.".dnp.annotated.maf";

    my $f_bam_list=$run_dir."/input.bam.list";

    open(OUTB,">$f_bam_list"); 

    foreach my $s (`ls $run_dir`) 
	{

	my $str=$s; 
	chomp($str);
 
	my $dir_2=$run_dir."/".$str; 

	if(-d $dir_2) 
	{
 	my $f_bam=$dir_2."/".$str.".remDup.bam";

	if(-e $f_bam) { print OUTB $str,"_T","\t",$f_bam,"\n"; }
  				
	}
			
	}
	
    open(DNP, ">$job_files_dir/$current_job_file") or die $!;
    print DNP "#!/bin/bash\n";
    ## remove snv nearby an indel ##
    print DNP "      ".$run_script_path."remove_nearby_snv.pl $f_maf $f_maf_rm_snv"."\n";
    ## annotate dnp ##
    print DNP "      ".$run_script_path."cocoon.pl $f_maf_rm_snv $f_maf_dnp_tmp $log_dir --bam $f_bam_list --samt $samtoolexe --merge --genome $h38_REF --gtf $f_gtf --snvonly"."\n";
    ## add dnp to the maf ##
    print DNP "		 ".$run_script_path."add_dnp.pl $f_maf_rm_snv $f_maf_dnp_tmp_merge $f_maf_dnp"."\n";
    ### remove tmp files ##
    print DNP "rm $f_maf_dnp_tmp_merge\n";
    print DNP "rm $f_maf_dnp_tmp\n";
    print DNP "rm $f_maf_rm_snv\n"; 
    print DNP "rm $f_maf_removed\n";
	
    close DNP;

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>10000] rusage[mem=10000]\" -M 10000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;

    system ($bsub_com);

}
	

exit;

sub bsub_fq2bam{

    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

    $current_job_file = "j0_fq2bam_".$sample_name.".sh";
    my $IN_fq1 = $sample_full_path."/".$sample_name.".R1.fq.gz";
    my $IN_fq2 = $sample_full_path."/".$sample_name.".R2.fq.gz";
    my $out_bam=$sample_full_path."/".$sample_name.".bam";
    my $out_sorted_bam=$sample_full_path."/".$sample_name.".sorted.bam";
    my $out_sorted_bam_bai=$sample_full_path."/".$sample_name.".sorted.bam.bai";
    my $out_rem_bam=$sample_full_path."/".$sample_name.".remDup.bam";
    my $out_metrics=$sample_full_path."/".$sample_name.".remdup.metrics.txt";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;

    open(FQ2BAM, ">$job_files_dir/$current_job_file") or die $!;
    print FQ2BAM "#!/bin/bash\n";
    print "export PATH=/miniconda/envs/align_dnaseq/bin:\${PATH}","\n"; 
    print FQ2BAM "              ".$run_py_script_path."align_dnaseq.py --out-prefix output --cpu 40 --flowcell HFMFWDSXY --index-sequencer CCAGTAGCGT-ATGTATTGGC --known-sites $f_known_site --lane 2 --library-preparation TWCE-HT191P1-S1H1A3Y3D1_1-lib1 --platform ILLUMINA --reference $GENOME $IN_fq1 $IN_fq2","\n";;
    close FQ2BAM;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

}

sub bsub_mutect2{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

	my $out_mutect2=$sample_full_path."/mutect2";	
	if(!(-d $out_mutect2)) { `mkdir $out_mutect2`; } 
    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
 		
    foreach my $chr (@chrlist)
    {
        my $chr1=$chr;

        if($chr_status==1)
        {
        $chr1="chr".$chr;
        }

        $current_job_file = "j1_mutect2_".$sample_name."_".$chr1.".sh";
        my $IN_bam_T = $sample_full_path."/".$sample_name.".remDup.bam";
        my $lsf_out=$lsf_file_dir."/".$current_job_file."_".$chr1.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file."_".$chr1.".err";
	my $f_out_gz=$out_mutect2."/".$chr1."-f1r2.tar.gz"; 
	my $f_out_vcf=$out_mutect2."/".$chr1."-unfiltered.vcf";

        `rm $lsf_out`;
        `rm $lsf_err`;

        open(MUTECT2, ">$job_files_dir/$current_job_file") or die $!;
        print MUTECT2 "#!/bin/bash\n";
		print MUTECT2 "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK Mutect2 -I $IN_bam_T -R $h38_REF -L $chr1 --germline-resource $GNOMAD_VCF -pon $PANEL_OF_NORMALS_VCF --f1r2-tar-gz $f_out_gz -O $f_out_vcf","\n";
		close MUTECT2;
    	my $sh_file=$job_files_dir."/".$current_job_file;

    	$bsub_com = "bsub -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);

    }

}

sub bsub_filter_mutect2 {

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j2_filter_mutect2_".$sample_name.".sh";
	my $out_mutect2=$sample_full_path."/mutect2";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
	my $all_unfiltered_input=""; 
    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    foreach my $chr (@chrlist)
    {
        my $chr1=$chr;

        if($chr_status==1)
        {
        $chr1="chr".$chr;
        }

	if($all_unfiltered_input eq "") { $all_unfiltered_input="I=".$out_mutect2."/".$chr1."-unfiltered.vcf"; } 
	else { $all_unfiltered_input=$all_unfiltered_input." "."I=".$out_mutect2."/".$chr1."-unfiltered.vcf";}
	}

	my $all_unfiltered_stats_input=""; 	
    foreach my $chr (@chrlist)
    {
        my $chr1=$chr;

        if($chr_status==1)
        {
        $chr1="chr".$chr;
        }

    	if($all_unfiltered_stats_input eq "") { $all_unfiltered_stats_input="-stats ".$out_mutect2."/".$chr1."-unfiltered.vcf.stats"; }
    	else { $all_unfiltered_stats_input=$all_unfiltered_stats_input." "."-stats ".$out_mutect2."/".$chr1."-unfiltered.vcf.stats";}
    }

	my $all_f1r2_input="";

    foreach my $chr (@chrlist)
    {
        my $chr1=$chr;

        if($chr_status==1)
        {
        $chr1="chr".$chr;
        }

        if($all_f1r2_input eq "") { $all_f1r2_input="-I ".$out_mutect2."/".$chr1."-f1r2.tar.gz"; }
        else { $all_f1r2_input=$all_f1r2_input." "."-I ".$out_mutect2."/".$chr1."-f1r2.tar.gz";}
    }

	my $IN_bam_T = $sample_full_path."/".$sample_name.".remDup.bam";
	my $f_merged_vcf=$out_mutect2."/merged-unfiltered.vcf";
    my $f_merged_status=$out_mutect2."/merged-unfiltered.vcf.stats"; 
	my $f_ori=$out_mutect2."/read-orientation-model.tar.gz";
	my $f_sum=$out_mutect2."/getpileupsummaries.table"; 
	my $f_tab=$out_mutect2."/contamination.table"; 
	my $f_tab2=$out_mutect2."/contamination2.table";

	my $f_seg=$out_mutect2."/segments.table";	
   	my $f_filtered_vcf=$out_mutect2."/filtered.vcf"; 
 
	open(MUTECT2F, ">$job_files_dir/$current_job_file") or die $!;

        print MUTECT2F "#!/bin/bash\n";
	print MUTECT2F "$java_bin -Xmx16g -jar $picardexe GatherVcfs $all_unfiltered_input O=$f_merged_vcf","\n";
	print MUTECT2F "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK MergeMutectStats $all_unfiltered_stats_input -O $f_merged_status\n"; 
	print MUTECT2F  "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK LearnReadOrientationModel $all_f1r2_input -O $f_ori\n";
        print MUTECT2F  "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx8g -jar $GATK GetPileupSummaries -I $IN_bam_T -V $COMMON_BIALLELIC -L $COMMON_BIALLELIC -O $f_sum\n";
	print MUTECT2F  "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK CalculateContamination -I $f_sum --tumor-segmentation $f_seg -O $f_tab\n";
	print MUTECT2F  "		".$run_script_path."check_contamination.pl $f_tab $f_tab2","\n"; 
	print MUTECT2F  "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK FilterMutectCalls -V $f_merged_vcf -R $h38_REF --tumor-segmentation $f_seg --contamination-table $f_tab2 --ob-priors $f_ori -O $f_filtered_vcf\n";	
	close MUTECT2F;

 	my $sh_file=$job_files_dir."/".$current_job_file;

        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    
        print $bsub_com;
        system ($bsub_com);

}

sub bsub_parse_mutect2{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j3_parse_mutect2_".$sample_name.".sh";

	my $out_mutect2=$sample_full_path."/mutect2";
    my $f_filtered_vcf=$out_mutect2."/filtered.vcf";
	my $f_passed_vcf=$out_mutect2."/filtered.pass.vcf";
	my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
	
    open(MUTECT2P, ">$job_files_dir/$current_job_file") or die $!;
    print MUTECT2P "#!/bin/bash\n";
    print MUTECT2P "export JAVA_HOME=$java_dir\n";
    print MUTECT2P "export JAVA_OPTS=\"-Xmx10g\"\n";
    print MUTECT2P "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print MUTECT2P "cat > $out_mutect2/mutect2_dbsnp_filter.input <<EOF\n";
    print MUTECT2P "mutect.dbsnp.annotator = $snpsift\n";
    print MUTECT2P "streka.dbsnp.db = $DB_SNP_NO_COSMIC\n";
    print MUTECT2P "streka.dbsnp.rawvcf = ./filtered.pass.vcf\n";
    print MUTECT2P "streka.dbsnp.mode = filter\n";
    print MUTECT2P "streka.dbsnp.passfile  = ./mutect2.somatic.dbsnp_pass.vcf\n";
    print MUTECT2P "streka.dbsnp.dbsnpfile = ./mutect2.somatic.dbsnp_present.vcf\n";
    print MUTECT2P "EOF\n";
    print MUTECT2P "cd $out_mutect2\n";
    print MUTECT2P "		".$run_script_path."filter_mutect2.pl $f_filtered_vcf $f_passed_vcf $mincov_t $minvaf\n";
    print MUTECT2P "     ".$run_script_path."dbsnp_filter.pl ./mutect2_dbsnp_filter.input\n";
    close MUTECT2P;	

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub  -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\'  -o $lsf_out -e $lsf_err bash $sh_file\n";
    
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


    $current_job_file = "j4_vcf_2_maf.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    #my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
	my $f_m_vcf_in=$sample_full_path."/mutect2/"."mutect2.somatic.dbsnp_pass.vcf";

    open(MAF, ">$job_files_dir/$current_job_file") or die $!;
    print MAF "#!/bin/bash\n";
    #print MAF "#BSUB -n 1\n";
    #print MAF "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print MAF "#BSUB -M 30000000\n";
    #print MAF "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print MAF "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print MAF "#BSUB -J $current_job_file\n";
    #print MAF "#BSUB -q ding-lab\n";
    #print MAF "#BSUB -a \'docker(scao/dailybox)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    #print MAF "#BSUB -q research-hpc\n";
    #print MAF "#BSUB -w \"$hold_job_file\"","\n";
  	
    print MAF "F_VCF_1=".$sample_full_path."/merged.mutect2.vcf\n";
    print MAF "F_VCF_2=".$sample_full_path."/".$sample_name.".mutect2.vcf\n";
    print MAF "F_VEP_1=".$sample_full_path."/merged.VEP.mutect2.vcf\n";
    print MAF "F_VEP_2=".$sample_full_path."/".$sample_name.".mutect2.vep.vcf\n";
    print MAF "F_maf=".$sample_full_path."/".$sample_name.".mutect2.maf\n";
    print MAF "RUNDIR=".$sample_full_path."\n";
	
    print MAF "F_log=".$sample_full_path."/vep.merged.mutect2.log"."\n";
    print MAF "cat > \${RUNDIR}/vep.merged.mutect2.input <<EOF\n";
    print MAF "merged.vep.vcf = ./merged.mutect2.vcf\n";
    print MAF "merged.vep.output = ./merged.VEP.mutect2.vcf\n";
    print MAF "merged.vep.vep_cmd = $vepannot\n";
    print MAF "merged.vep.cachedir = $vepcache\n";
    print MAF "merged.vep.reffasta = $f_ref_annot\n";
    print MAF "merged.vep.assembly = GRCh38\n";
    print MAF "EOF\n";
    print MAF "rm \${F_log}\n";
	
  
	### vep and vcf2maf annotation for all variants to get the annotated gene name for each variant ##
    #print MAF "export PERL5LIB=/storage1/fs1/dinglab/Active/Projects/austins2/test/perl-5.24.1/lib/perl5:/storage1/fs1/dinglab/Active/Projects/austins2/test/perl-5.24.1/x86_64-linux-gnu-thread-multi/:\$PERL5LIB\n";
    #print MAF ". $script_dir/set_envvars\n";
    print MAF "cd $sample_full_path\n";
    #print MAF ". $script_dir/set_envvars\n";
    print MAF "cp $f_m_vcf_in \${F_VCF_1}","\n";  
    print MAF "     ".$run_script_path."vep_annotator.pl ./vep.merged.mutect2.input >&./vep.merged.mutect2.log\n";
    print MAF "rm \${F_VCF_2}\n";
    print MAF "rm \${F_VEP_2}\n";
    print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
    print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
    print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot --file-tsl $TSL_DB\n";
	
    close MAF;

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'  -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

   # $bsub_com = "bsub -G $compute_group -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
 
   #print $bsub_com;
#	system ($bsub_com);

 }

 
