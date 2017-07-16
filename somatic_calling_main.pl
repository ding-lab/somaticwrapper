######### Song Cao###########
## pipeline for somatic variant callings ##
#	somatic_variant_callings.pl #
###	updated date: 04/05/2017 ###
### updated date: 04/18/2017 ###
### add vcf2maf.pl ###
### 07/14/2017 ##
### add mutect and cosmic databse for improving filtering ##

#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
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
This script will process rna-seq data for TCGA samples. 
Pipeline version: $version
$yellow     Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$green       [1]  Run streka
$red 		 [2]  Run Varscan
$normal 	 [3]  Run Mutect
$yellow 	 [4]  Parse streka result
$purple 	 [5]  Parse VarScan result
$normal      [6]  Parse Mutect result
$cyan 		 [7]  Run Pindel
$gray 		 [8]  Run VEP annotation
$gray 		 [9]  Parse Pindel
$gray 		 [10]  Merge vcf files  
$gray 		 [11] generate maf file 
$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}
die $usage unless ($step_number >=0)&&(($step_number <= 10));
my $email = "scao\@wustl\.edu";
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];
my $HOME1="/gscmnt/gc2524/dinglab";
#store job files here
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
my $script_dir="/gscuser/scao/scripts/git/somaticwrapper";
#my $run_script_path = `dirname $0`;
my $run_script_path=$script_dir; 
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";
print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
my $mutect="/gscuser/rmashl/Software/bin/gatk/3.7/GenomeAnalysisTK.jar";
my $STRELKA_DIR="/gscmnt/gc2525/dinglab/rmashl/Software/bin/strelka/1.0.14/bin";
my $h37_REF="/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37V1/GRCh37-lite-chr_with_chrM.fa";
my $f_exac="/gscmnt/gc2741/ding/qgao/tools/vcf2maf-1.6.11/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
my $h37_REF_bai="/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37/GRCh37-lite-chr_with_chrM.fa.fai";
my $pindel="/gscuser/qgao/tools/pindel/pindel";
my $PINDEL_DIR="/gscuser/qgao/tools/pindel";
my $gatk="/gscuser/scao/tools/GenomeAnalysisTK.jar";
my $f_centromere="/gscmnt/gc3015/dinglab/medseq/Jiayin_Germline_Project/PCGP/data/pindel-centromere-exclude.bed";

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure
#&check_input_dir($run_dir);
# start data processsing

if ($step_number < 10) {
    #begin to process each sample
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
                $current_job_file="";
                if($step_number==0)
                {  
				   &bsub_strelka();
				   &bsub_varscan();
				   &bsub_mutect();
				   &bsub_parse_strelka();
				   &bsub_parse_varscan();
				   &bsub_pindel();
				   &bsub_vep();
				   &bsub_parse_pindel();
				   &bsub_merge_vcf();
				   &bsub_vcf_2_maf();
				  # &bsub_pindel();	
				}
                 elsif ($step_number == 1) {
                    &bsub_strelka();
                } elsif ($step_number == 2) {
                    &bsub_varscan(1);
                } elsif ($step_number == 3) {
                    &bsub_mutect(1);
                }
				elsif ($step_number == 4) {
                    &bsub_parse_strelka(1);
                } 
				elsif ($step_number == 5) {
                    &bsub_parse_varscan(1);
                }elsif ($step_number == 6) {
                    &bsub_pindel(1);
                }elsif ($step_number == 7) {
                    &bsub_vep(1);
                }elsif ($step_number == 8) {
                    &bsub_parse_pindel(1);
                }elsif ($step_number == 9) {
                    &bsub_merge_vcf(1);
                }elsif ($step_number == 10) {
                    &bsub_vcf_2_maf(1);
                }

				# how to run pindel #
				# /gscmnt/gc8001/info/model_data/8ad63791a648435b95a4be23ac3a18a5/buildb8e7ade8c3244bd483c274edeceeb8b7/alignments/cfa9c81317e346f5aa458e5e9220f41e.bam 500 TCGA-EX-A1H5-01A-31D-A34H-09 #
#/gscmnt/gc8001/info/model_data/c777ae8a6fdb49b7bf400056d15a6788/build8057e74cc9f64a50a13c33763af10116/alignments/b127fcaaecc54abb833316fa3153e492.bam 500 TCGA-EX-A1H5-10A-01D-A200-09#	
# bsub -q ding-lab -M 16000000 -n 4 -R 'span[hosts=1] rusage[mem=16000]' -oo /gscmnt/gc2532/dinglab/scao/pindel/CESC/pindel-logs/TCGA-C5-A0TN.chr9.log /gscuser/kye/gc2532/projects/PCGP_rerun/pindel -T 4 -f /gscmnt/gc2532/dinglab/projects/PCGP_rerun/GRCh37-lite.fa -i /gscmnt/gc2532/dinglab/scao/pindel/CESC/pindel-configs/TCGA-C5-A0TN.config -o /gscmnt/gc2532/dinglab/scao/pindel/CESC/pindel-outputs/TCGA-C5-A0TN/TCGA-C5-A0TN.chr9 -c 9 -m 6 -w 1 -J /gscmnt/gc3015/dinglab/medseq/Jiayin_Germline_Project/PCGP/data/pindel-centromere-exclude.bed #
 
           }
        }
    }
}

#######################################################################
# send email to notify the finish of the analysis
if (($step_number == 0) || ($step_number == 10)) {
    print $yellow, "Submitting the job for sending an email when the run finishes ",$sample_name, "...",$normal, "\n";
    $hold_job_file = $current_job_file;
    $current_job_file = "Email_run_".$$.".sh";
    open(EMAIL, ">$job_files_dir/$current_job_file") or die $!;
    print EMAIL "#!/bin/bash\n";
    print EMAIL "#BSUB -n 1\n";
    print EMAIL "#BSUB -o $lsf_file_dir","\n";
    print EMAIL "#BSUB -e $lsf_file_dir","\n";
    print EMAIL "#BSUB -J $current_job_file\n";
    print EMAIL "#BSUB -w \"$hold_job_file\"","\n";
    print EMAIL $run_script_path."send_email.pl ".$run_dir." ".$email."\n";
    close EMAIL;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}
#######################################################################
if ($step_number == 0) {
    print $green, "All jobs are submitted! You will get email notification when this run is completed.\n",$normal;
}

exit;


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
    open(STREKA, ">$job_files_dir/$current_job_file") or die $!;
    print STREKA "#!/bin/bash\n";
    print STREKA "#BSUB -n 1\n";
    print STREKA "#BSUB -R \"rusage[mem=30000]\"","\n";
    print STREKA "#BSUB -M 30000000\n";
    print STREKA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print STREKA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print STREKA "#BSUB -J $current_job_file\n";
	print STREKA "scr_t0=\`date \+\%s\`\n";
    print STREKA "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print STREKA "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print STREKA "myRUNDIR=".$sample_full_path."/strelka\n";
	print STREKA "STATUSDIR=".$sample_full_path."/status\n";
    print STREKA "RESULTSDIR=".$sample_full_path."/results\n";
	print STREKA "SG_DIR=".$sample_full_path."/strelka\n"; 
	print STREKA "RUNDIR=".$sample_full_path."\n";
	print STREKA "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
   	print STREKA "CONFDIR="."/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/config\n";
 	print STREKA "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
	print STREKA "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
	print STREKA "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
	print STREKA "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print STREKA "if [ ! -d \${myRUNDIR} ]\n";
	print STREKA "then\n";
	print STREKA "mkdir \${myRUNDIR}\n";
	print STREKA "fi\n";
    print STREKA "if [ -d \${STRELKA_OUT} ]\n";
    print STREKA "then\n";
    print STREKA "rm -rf \${STRELKA_OUT}\n";
    print STREKA "fi\n";
	print STREKA "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n"; 
   	print STREKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
	print STREKA "else\n";
   	print STREKA "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
	print STREKA "fi\n";
	print STREKA "put_cmd=\"ln -s\"\n";
	print STREKA "del_cmd=\"rm -f\"\n";
	print STREKA "del_local=\"rm -f\"\n";
	print STREKA "statfile=incomplete.strelka\n";
	print STREKA "localstatus=".$sample_full_path."/status/\$statfile\n";
	print STREKA "touch \$localstatus\n";
	print STREKA "   ".$STRELKA_DIR."/configureStrelkaWorkflow.pl --normal \$NBAM --tumor \$TBAM --ref ". $h37_REF." --config \$CONFDIR/strelka.ini --output-dir \$STRELKA_OUT\n";
	print STREKA "cd \$STRELKA_OUT\n";
	print STREKA "make -j 16\n";
    close STREKA;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}

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
    print VARSCAN "#BSUB -n 1\n";
    print VARSCAN "#BSUB -R \"rusage[mem=30000]\"","\n";
    print VARSCAN "#BSUB -M 30000000\n";
    print VARSCAN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print VARSCAN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print VARSCAN "#BSUB -J $current_job_file\n";
	print VARSCAN "#BSUB -w \"$hold_job_file\"","\n";
  	print VARSCAN "scr_t0=\`date \+\%s\`\n";
    #print VARSCAN "chralt=\${chr\/:\/_}\n";
	#print VARSCAN "dir=\$chralt\n";
	print VARSCAN "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print VARSCAN "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print VARSCAN "myRUNDIR=".$sample_full_path."/varscan\n";
    print VARSCAN "STATUSDIR=".$sample_full_path."/status\n";
    print VARSCAN "RESULTSDIR=".$sample_full_path."/varscan_results\n";
    print VARSCAN "RUNDIR=".$sample_full_path."\n";
    #print VARSCAN "numgps=10\n";
	#print VARSCAN "SEQS=\"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y\"\n";
    print VARSCAN "CONFDIR="."/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/config\n";
   	print VARSCAN "GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip\n";
	print VARSCAN "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
	print VARSCAN "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print VARSCAN "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print VARSCAN "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
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
    print VARSCAN "statfile=incomplete.vs_som_snvindels\n";
	print VARSCAN "localstatus=\${RUNDIR}\/status\/\${statfile}\n";
	print VARSCAN "if [ ! -d \${myRUNDIR}\/status ]\n";
    print VARSCAN "then\n";
    print VARSCAN "mkdir \${myRUNDIR}\/status\n";
    print VARSCAN "fi\n";
	print VARSCAN "touch \${localstatus}\n";
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
	print VARSCAN "ncols=\$(echo \"3*( \$(wc -l < \$BAMLIST) +1)\"|bc)\n";
	print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somatic - \${TMPBASE} --mpileup 1 --p-value 0.99 --somatic-p-value 0.05 --min-coverage-normal 30 --min-coverage-tumor 22 --min-var-freq 0.08 --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --strand-filter 1 --min-avg-qual 15 --output-vcf 1 --output-snp \${snvoutbase} --output-indel \${indeloutbase} &> \${LOG}\n";
  	close VARSCAN;	
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}

sub bsub_mutect{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j3_mutect_".$sample_name.".sh";
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

    open(MUTECT, ">$job_files_dir/$current_job_file") or die $!;
    print MUTECT "#!/bin/bash\n";
    print MUTECT "#BSUB -n 1\n";
    print MUTECT "#BSUB -R \"rusage[mem=30000]\"","\n";
    print MUTECT "#BSUB -M 30000000\n";
    print MUTECT "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print MUTECT "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print MUTECT "#BSUB -J $current_job_file\n";
    print MUTECT "scr_t0=\`date \+\%s\`\n";
    print MUTECT "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print MUTECT "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print MUTECT "myRUNDIR=".$sample_full_path."/mutect\n";
	print MUTECT "rawvcf=".$sample_full_path."/mutect/mutect.raw.vcf\n";
	print MUTECT "rawvcfgvip=".$sample_full_path."/mutect/mutect.raw.gvip.vcf\n";
	print MUTECT "rawvcfsnv=".$sample_full_path."/mutect/mutect.raw.gvip.snv.vcf\n";
	print MUTECT "rawvcfindel=".$sample_full_path."/mutect/mutect.raw.gvip.indel.vcf\n";
    print MUTECT "RUNDIR=".$sample_full_path."\n";
    print MUTECT "CONFDIR="."/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/config\n";
    print MUTECT "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print MUTECT "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print MUTECT "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print MUTECT "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print MUTECT "if [ ! -d \${myRUNDIR} ]\n";
    print MUTECT "then\n";
    print MUTECT "mkdir \${myRUNDIR}\n";
    print MUTECT "fi\n";
    print MUTECT "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
    print MUTECT "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print MUTECT "else\n";
    print MUTECT "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print MUTECT "fi\n";
    print MUTECT "java  \${JAVA_OPTS} -jar $mutect  -R $h37_REF  -T MuTect2 -I:tumor $IN_bam_N -I:normal $IN_bam_T  -mbq  10  -rf DuplicateRead    -rf UnmappedRead    -stand_call_conf  10.0    -o  \${raw.vcf}\n";
    print MUTECT "     ".$run_script_path."genomevip_label.pl mutect \${rawvcf} \${rawvcfgvip}\n";
    print MUTECT "java \${JAVA_OPTS} -jar $mutect  -R $h37_REF  -T SelectVariants  -V  \${rawvcfgvip}  -o  \${rawvcfsnv}   -selectType SNP -selectType MNP\n";
	print MUTECT "java \${JAVA_OPTS} -jar $mutect  -R $h37_REF  -T SelectVariants  -V  \${rawvcfgvip} -o \${rawvcfindel}  -selectType INDEL\n";
	close MUTECT;
	$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}

sub bsub_parse_strelka{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j4_parse_strelka".$sample_name.".sh";

    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    open(STREKAP, ">$job_files_dir/$current_job_file") or die $!;

    print STREKAP "#!/bin/bash\n";
    print STREKAP "#BSUB -n 1\n";
    print STREKAP "#BSUB -R \"rusage[mem=30000]\"","\n";
    print STREKAP "#BSUB -M 30000000\n";
    print STREKAP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print STREKAP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print STREKAP "#BSUB -J $current_job_file\n";
	print STREKAP "#BSUB -w \"$hold_job_file\"","\n";
    print STREKAP "#BSUB -q ding-lab\n";
    print STREKAP "scr_t0=\`date \+\%s\`\n";
    print STREKAP "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print STREKAP "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print STREKAP "myRUNDIR=".$sample_full_path."/strelka\n";
    print STREKAP "STATUSDIR=".$sample_full_path."/status\n";
    print STREKAP "RESULTSDIR=".$sample_full_path."/results\n";
    print STREKAP "SG_DIR=".$sample_full_path."/strelka\n";
    print STREKAP "RUNDIR=".$sample_full_path."\n";
    print STREKAP "STRELKA_OUT=".$sample_full_path."/strelka/strelka_out"."\n";
	print STREKAP "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print STREKAP "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print STREKAP "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print STREKAP "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print STREKAP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print STREKAP "cat > \${myRUNDIR}/strelka_out/results/strelka_dbsnp_filter.snv.input <<EOF\n";
    print STREKAP "streka.dbsnp.snv.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar\n";
    print STREKAP "streka.dbsnp.snv.db = /gscmnt/gc3027/dinglab/medseq/cosmic/00-All.brief.pass.cosmic.vcf\n";
    print STREKAP "streka.dbsnp.snv.rawvcf = ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
    print STREKAP "streka.dbsnp.snv.mode = filter\n";
    print STREKAP "streka.dbsnp.snv.passfile  = ./strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
    print STREKAP "streka.dbsnp.snv.dbsnpfile = ./strelka.somatic.snv.all.gvip.dbsnp_present.vcf\n";
    print STREKAP "EOF\n";
  #	print STREKAP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.snv.input <<EOF\n";
  #  print STREKAP "strelka.vep.vcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
  #  print STREKAP "strelka.vep.output = ./varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf\n";
  #  print STREKAP "strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
  #  print STREKAP "strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
  #  print STREKAP "strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
 #   print STREKAP "strelka.vep.assembly = GRCh37\n";
 #   print STREKAP "EOF\n";
 #   print STREKAP "cat > \${RUNDIR}/varscan/vs_vep.indel.input <<EOF\n";
 #   print STREKAP "varscan.vep.vcf = ./varscan.out.som_indel.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
 #   print STREKAP "varscan.vep.output = ./varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf\n";
 #   print STREKAP "varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
 #   print STREKAP "varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
 #   print STREKAP "varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
  #  print STREKAP "varscan.vep.assembly = GRCh37\n";
  #  print STREKAP "EOF\n";  
	print STREKAP "cat > \${myRUNDIR}/strelka_out/results/strelka_dbsnp_filter.indel.input <<EOF\n";
   	print STREKAP "streka.dbsnp.indel.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar\n";
    print STREKAP "streka.dbsnp.indel.db = /gscmnt/gc3027/dinglab/medseq/cosmic/00-All.brief.pass.cosmic.vcf\n";
    print STREKAP "streka.dbsnp.indel.rawvcf = ./strelka.somatic.indel.strlk_pass.gvip.vcf\n";
    print STREKAP "streka.dbsnp.indel.mode = filter\n";
    print STREKAP "streka.dbsnp.indel.passfile  = ./strelka.somatic.indel.all.gvip.dbsnp_pass.vcf\n";
    print STREKAP "streka.dbsnp.indel.dbsnpfile = ./strelka.somatic.indel.all.gvip.dbsnp_present.vcf\n";
	print STREKAP "EOF\n";
	print STREKAP "FP_BAM=\`awk \'{if(NR==1){print \$1}}\' \${RUNDIR}/varscan/bamfilelist.inp\`\n";
	print STREKAP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_fpfilter.snv.input <<EOF\n";
	print STREKAP "strelka.fpfilter.snv.bam_readcount = /gscmnt/gc2525/dinglab/rmashl/Software/bin/bam-readcount/0.7.4/bam-readcount\n";
	print STREKAP "strelka.fpfilter.snv.bam_file = \${FP_BAM}\n";
	print STREKAP "strelka.fpfilter.snv.REF = $h37_REF\n"; 
	print STREKAP "strelka.fpfilter.snv.variants_file = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
	print STREKAP "strelka.fpfilter.snv.passfile = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_pass.vcf\n";
	print STREKAP "strelka.fpfilter.snv.failfile = \${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_fail.vcf\n";
	print STREKAP "strelka.fpfilter.snv.min_mapping_qual = 0\n";
	print STREKAP "strelka.fpfilter.snv.min_base_qual = 15\n";
	print STREKAP "strelka.fpfilter.snv.min_num_var_supporting_reads = 4\n";
	print STREKAP "strelka.fpfilter.snv.min_var_allele_freq = 0.05\n";
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
	print STREKAP "     ".$run_script_path."genomevip_label.pl Strelka ./all.somatic.snvs.vcf ./strelka.somatic.snv.all.gvip.vcf\n";
    print STREKAP "     ".$run_script_path."genomevip_label.pl Strelka ./all.somatic.indels.vcf ./strelka.somatic.indel.all.gvip.vcf\n";
	print STREKAP " 	"."cp ./strelka.somatic.snv.all.gvip.vcf ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
	print STREKAP "     "."cp ./strelka.somatic.indel.all.gvip.vcf ./strelka.somatic.indel.strlk_pass.gvip.vcf\n";
   # print STREKAP "     ".$run_script_path."genomevip_label.pl Strelka ./passed.somatic.snvs.vcf ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
   # print STREKAP "     ".$run_script_path."genomevip_label.pl Strelka ./passed.somatic.indels.vcf ./strelka.somatic.indel.strlk_pass.gvip.vcf\n"; 
	print STREKAP "     ".$run_script_path."dbsnp_filter.pl ./strelka_dbsnp_filter.snv.input\n";
    print STREKAP "     ".$run_script_path."dbsnp_filter.pl ./strelka_dbsnp_filter.indel.input\n";
    print STREKAP "     ".$run_script_path."snv_filter.pl ./strelka_fpfilter.snv.input\n";  
    #print STREKAP "     ".$run_script_path."vep_annotator.pl ./vs_vep.snv.input >& ./vs_vep.snv.log\n";
    #print STREKAP "     ".$run_script_path."vep_annotator.pl ./vs_vep.indel.input >& ./vs_vep.indel.log\n";	
	close STREKAP;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com ); 
# print STREKAP "EOF\n";

}
sub bsub_parse_varscan{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


  	$current_job_file = "j5_parse_varscan".$sample_name.".sh";

    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    open(VARSCANP, ">$job_files_dir/$current_job_file") or die $!;

    print VARSCANP "#!/bin/bash\n";
    print VARSCANP "#BSUB -n 1\n";
    print VARSCANP "#BSUB -R \"rusage[mem=30000]\"","\n";
    print VARSCANP "#BSUB -M 30000000\n";
    print VARSCANP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print VARSCANP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print VARSCANP "#BSUB -J $current_job_file\n"; 
    print VARSCANP "#BSUB -q ding-lab\n";
	print VARSCANP "#BSUB -w \"$hold_job_file\"","\n";
    print VARSCANP "scr_t0=\`date \+\%s\`\n";
	print VARSCANP "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print VARSCANP "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print VARSCANP "myRUNDIR=".$sample_full_path."/varscan\n";
    print VARSCANP "STATUSDIR=".$sample_full_path."/status\n";
    print VARSCANP "RUNDIR=".$sample_full_path."\n";
    print VARSCANP "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print VARSCANP "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print VARSCANP "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print VARSCANP "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print VARSCANP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print VARSCANP "cat > \${RUNDIR}/varscan/vs_dbsnp_filter.snv.input <<EOF\n";
	print VARSCANP "varscan.dbsnp.snv.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar\n";
	print VARSCANP "varscan.dbsnp.snv.db = /gscmnt/gc3027/dinglab/medseq/cosmic/00-All.brief.pass.cosmic.vcf\n";
	print VARSCANP "varscan.dbsnp.snv.rawvcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.snv.mode = filter\n";
	print VARSCANP "varscan.dbsnp.snv.passfile  = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.snv.dbsnpfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf\n";
	print VARSCANP "EOF\n";
 # 	print VARSCANP "cat > \${RUNDIR}/varscan/vs_cosmic_check.snv.input <<EOF\n";
 #  print VARSCANP "varscan.dbsnp.snv.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar\n";
 #  print VARSCANP "varscan.dbsnp.snv.db = /gscmnt/gc2525/dinglab/rmashl/Software/bin/dbSNP/NCBI/snp142/GRCh37/00-All.brief.vcf\n";
 #  print VARSCANP "varscan.dbsnp.snv.rawvcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf\n";
 #  print VARSCANP "varscan.dbsnp.snv.mode = filter\n";
 #  print VARSCANP "varscan.dbsnp.snv.passfile  = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
 #  print VARSCANP "varscan.dbsnp.snv.dbsnpfile = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf\n";
 #  print VARSCANP "EOF\n";
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_dbsnp_filter.indel.input <<EOF\n";
	print VARSCANP "varscan.dbsnp.indel.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar\n";
	print VARSCANP "varscan.dbsnp.indel.db = /gscmnt/gc3027/dinglab/medseq/cosmic/00-All.brief.pass.cosmic.vcf\n";
	print VARSCANP "varscan.dbsnp.indel.rawvcf = ./varscan.out.som_indel.gvip.Somatic.hc.vcf\n";
	print VARSCANP "varscan.dbsnp.indel.mode = filter\n";
	print VARSCANP "varscan.dbsnp.indel.passfile  = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.dbsnp.indel.dbsnpfile = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_present.vcf\n";
	print VARSCANP "EOF\n";
	print VARSCANP "FP_BAM=\`awk \'{if(NR==2){print \$1}}\' \${RUNDIR}/varscan/bamfilelist.inp\`\n";	
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_fpfilter.somatic.snv.input <<EOF\n";
	print VARSCANP "varscan.fpfilter.snv.bam_readcount = /gscmnt/gc2525/dinglab/rmashl/Software/bin/bam-readcount/0.7.4/bam-readcount\n";
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
	print VARSCANP "cat > \${RUNDIR}/varscan/vs_vep.snv.input <<EOF\n";
	print VARSCANP "varscan.vep.vcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print VARSCANP "varscan.vep.output = ./varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf\n";
	print VARSCANP "varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
	print VARSCANP "varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
	print VARSCANP "varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
	print VARSCANP "varscan.vep.assembly = GRCh37\n";
	print VARSCANP "EOF\n";
    print VARSCANP "cat > \${RUNDIR}/varscan/vs_vep.indel.input <<EOF\n";
    print VARSCANP "varscan.vep.vcf = ./varscan.out.som_indel.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
    print VARSCANP "varscan.vep.output = ./varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf\n";
    print VARSCANP "varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VARSCANP "varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VARSCANP "varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VARSCANP "varscan.vep.assembly = GRCh37\n";
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
    print VARSCANP "statfile=incomplete.vs_som_snvindels\n";
    print VARSCANP "localstatus=\${RUNDIR}\/status\/\${statfile}\n";
    print VARSCANP "if [ ! -d \${myRUNDIR}\/status ]\n";
    print VARSCANP "then\n";
    print VARSCANP "mkdir \${myRUNDIR}\/status\n";
    print VARSCANP "fi\n";
    print VARSCANP "touch \${localstatus}\n";
    print VARSCANP "cd \${RUNDIR}/varscan\n";
    print VARSCANP "TMPBASE=.\/varscan.out.som\n";
    print VARSCANP "LOG=\${TMPBASE}.log\n";
    print VARSCANP "snvoutbase=\${TMPBASE}_snv\n";
    print VARSCANP "indeloutbase=\${TMPBASE}_indel\n";
	print VARSCANP "cd \${RUNDIR}/varscan\n";
	print VARSCANP "     ".$run_script_path."genomevip_label.pl VarScan \${snvoutbase}.vcf  \${snvoutbase}.gvip.vcf\n";
	print VARSCANP "     ".$run_script_path."genomevip_label.pl VarScan \${indeloutbase}.vcf \${indeloutbase}.gvip.vcf\n";
	print VARSCANP "echo \'APPLYING PROCESS FILTER TO SOMATIC SNVS:\' &>> \${LOG}\n";
	print VARSCANP "mysnvorig=./\${snvoutbase}.gvip.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar processSomatic \${mysnvorig} --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_path."extract_somatic_other.pl <  \${mysnvorig}  > \${mysnvorig/%vcf/other.vcf}\n";
    print VARSCANP "for kk in Somatic Germline LOH ; do\n";
   	print VARSCANP "thisorig=\${mysnvorig/%vcf/\$kk.vcf}\n";
    print VARSCANP "thispass=\${mysnvorig/%vcf/\$kk.hc.vcf}\n";
   	print VARSCANP "thisfail=\${mysnvorig/%vcf/\$kk.lc.vcf}\n";
    print VARSCANP "     ".$script_dir."/extract_fail.sh  ./\${thisorig}  ./\${thispass}  ./\${thisfail}\n";
   	print VARSCANP "     ".$script_dir."/set_vcf_filter_label.sh  ./\${thisfail}  hc_fail\n";
	print VARSCANP "done\n";
	print VARSCANP "echo \'APPLYING PROCESS FILTER TO SOMATIC INDELS:\' &>> \$LOG\n";
	print VARSCANP "myindelorig=./\$indeloutbase.gvip.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar processSomatic \${myindelorig}   --min-tumor-freq  0.10   --max-normal-freq  0.05   --p-value  0.07  &>> \${LOG}\n";
	print VARSCANP "     ".$run_script_path."extract_somatic_other.pl <  \${myindelorig}  >  \${myindelorig/%vcf/other.vcf}\n";
	print VARSCANP "for kk in Somatic Germline LOH ; do\n";
    print VARSCANP "thisorig=\${myindelorig/%vcf/\$kk.vcf}\n";
   	print VARSCANP "thispass=\${myindelorig/%vcf/\$kk.hc.vcf}\n";
   	print VARSCANP "thisfail=\${myindelorig/%vcf/\$kk.lc.vcf}\n";
   	print VARSCANP "     ".$script_dir."/extract_fail.sh  ./\${thisorig}  ./\${thispass}  ./\${thisfail}\n";
   	print VARSCANP "     ".$script_dir."/set_vcf_filter_label.sh  ./\${thisfail}  hc_fail\n";
	print VARSCANP "done\n";
	print VARSCANP "scr_tf=\`date +%s\`\n";
	print VARSCANP "scr_dt=\$((scr_tf - scr_t0))\n";
	print VARSCANP "echo GVIP_TIMING_VARSCAN_DISCOVERY=\${scr_t0},\${scr_dt}\n";
	print VARSCANP "scr_t0=\${scr_tf}\n";
	print VARSCANP "echo \'APPLYING SOMATIC FILTER:\' &>> \${LOG}\n";
	print VARSCANP "thissnvorig=\${snvoutbase}.gvip.Somatic.hc.vcf\n";
	print VARSCANP "myindelorig=\${indeloutbase}.gvip.vcf\n";
	print VARSCANP "thissnvpass=\${snvoutbase}.gvip.Somatic.hc.somfilter_pass.vcf\n";
	print VARSCANP "thissnvfail=\${snvoutbase}.gvip.Somatic.hc.somfilter_fail.vcf\n";
	print VARSCANP "java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar somaticFilter  ./\${thissnvorig} --min-coverage  30   --min-reads2  4   --min-strands2  1   --min-avg-qual  20   --min-var-freq  0.10   --p-value  0.05   --indel-file  ./\${myindelorig} --output-file  ./\${thissnvpass}  &>> \${LOG}\n";
	print VARSCANP "     ".$script_dir."/extract_fail.sh ./\${thissnvorig}  ./\${thissnvpass}   ./\${thissnvfail}\n";
	print VARSCANP "     ".$script_dir."/set_vcf_filter_label.sh  ./\${thissnvfail}   somfilter_fail\n";
# $del_local  ./$thissnvorig
# $del_local  ./$myindelorig
	print VARSCANP "     ".$run_script_path."dbsnp_filter.pl  \${RUNDIR}/varscan/vs_dbsnp_filter.snv.input\n";
# $del_local ./varscan.out.som_snv.group0.chr1.gvip.Somatic.hc.somfilter_pass.vcf
	print VARSCANP "     ".$run_script_path."dbsnp_filter.pl \${RUNDIR}/varscan/vs_dbsnp_filter.indel.input\n";
# $del_local ./varscan.out.som_indel.group0.chr1.gvip.Somatic.hc.vcf
	print VARSCANP "     ".$run_script_path."snv_filter.pl  \${RUNDIR}/varscan/vs_fpfilter.somatic.snv.input\n";
	print VARSCANP "     ".$run_script_path."vep_annotator.pl ./vs_vep.snv.input >& ./vs_vep.snv.log\n";
	print VARSCANP "     ".$run_script_path."vep_annotator.pl ./vs_vep.indel.input >& ./vs_vep.indel.log\n";
# $del_local  ./varscan.out.som_snv.group0.chr1.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf
#mkdir -p $myRESULTSDIR
	print VARSCANP "scr_tf=\`date +%s\`\n";
	print VARSCANP "scr_dt=\$((scr_tf - scr_t0))\n";
	print VARSCANP "echo GVIP_TIMING_VARSCAN_FILTERING=\${scr_t0},\${scr_dt}\n";
	close VARSCANP; 
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

	}

sub bsub_pindel{
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

	$current_job_file = "j6_pindel".$sample_name.".sh";  
	my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    open(PINDEL, ">$job_files_dir/$current_job_file") or die $!;
    print PINDEL "#!/bin/bash\n";
    print PINDEL "#BSUB -n 4\n";
    print PINDEL "#BSUB -R \"span[hosts=1] rusage[mem=30000]\"","\n";
    print PINDEL "#BSUB -M 30000000\n";
    print PINDEL "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print PINDEL "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print PINDEL "#BSUB -J $current_job_file\n";
    print PINDEL "#BSUB -w \"$hold_job_file\"","\n";
    print PINDEL "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print PINDEL "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
	print PINDEL "myRUNDIR=".$sample_full_path."/pindel\n";
	print PINDEL "CONFIG=\${myRUNDIR}"."/".$sample_name.".config\n";
	print PINDEL "if [ ! -d \${myRUNDIR} ]\n";
    print PINDEL "then\n";
    print PINDEL "mkdir \${myRUNDIR}\n";
    print PINDEL "fi\n";
	print PINDEL "echo \"$IN_bam_T\t500\t$sample_name.T\" > \${CONFIG}\n";
    print PINDEL "echo \"$IN_bam_N\t500\t$sample_name.N\" >> \${CONFIG}\n";
	print PINDEL "$pindel -T 4 -f $h37_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"." -m 6 -w 1 -J $f_centromere\n";
	close PINDEL;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );	
	}

sub bsub_vep{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j7_vep".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    open(VEP, ">$job_files_dir/$current_job_file") or die $!; 
	print VEP "#!/bin/bash\n";
    print VEP "#BSUB -n 1\n";
    print VEP "#BSUB -R \"rusage[mem=30000]\"","\n";
    print VEP "#BSUB -M 30000000\n";
    print VEP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print VEP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print VEP "#BSUB -J $current_job_file\n";
    print VEP "#BSUB -q ding-lab\n";
    print VEP "#BSUB -w \"$hold_job_file\"","\n";
    print VEP "scr_t0=\`date \+\%s\`\n";
    print VEP "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print VEP "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print VEP "myRUNDIR=".$sample_full_path."/varscan\n";
    print VEP "STATUSDIR=".$sample_full_path."/status\n";
    print VEP "RUNDIR=".$sample_full_path."\n";
    print VEP "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print VEP "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print VEP "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print VEP "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print VEP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print VEP "cat > \${RUNDIR}/varscan/vs_vep.snv.input <<EOF\n";
    print VEP "varscan.vep.vcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
    print VEP "varscan.vep.output = ./varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf\n";
    print VEP "varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VEP "varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VEP "varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VEP "varscan.vep.assembly = GRCh37\n";
    print VEP "EOF\n";
    print VEP "cat > \${RUNDIR}/varscan/vs_vep.indel.input <<EOF\n";
    print VEP "varscan.vep.vcf = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf\n";
    print VEP "varscan.vep.output = ./varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf\n";
    print VEP "varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VEP "varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VEP "varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VEP "varscan.vep.assembly = GRCh37\n";
    print VEP "EOF\n";
  	print VEP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.snv.input <<EOF\n";
  	print VEP "strelka.vep.vcf = ./strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
    print VEP "strelka.vep.output = ./strelka.somatic_snv.current_final.gvip.Somatic.VEP.vcf\n";
  	print VEP "strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
  	print VEP "strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
  	print VEP "strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
 	print VEP "strelka.vep.assembly = GRCh37\n";
 	print VEP "EOF\n";
    print VEP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.indel.input <<EOF\n";
    print VEP "strelka.vep.vcf = ./strelka.somatic.indel.all.gvip.dbsnp_pass.vcf\n";
    print VEP "strelka.vep.output = ./strelka.somatic_indel.current_final.gvip.Somatic.VEP.vcf\n";
    print VEP "strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VEP "strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VEP "strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VEP "strelka.vep.assembly = GRCh37\n";
    print VEP "EOF\n";
	### VEP annotation for the initial vcf files ###
    print VEP "cat > \${RUNDIR}/varscan/vs_vep.snv.inital.input <<EOF\n";
    print VEP "varscan.vep.vcf = ./varscan.out.som_snv.gvip.vcf\n";
    print VEP "varscan.vep.output = ./varscan.out.som_snv.gvip.VEP.vcf\n";
    print VEP "varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VEP "varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VEP "varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VEP "varscan.vep.assembly = GRCh37\n";
    print VEP "EOF\n";
    print VEP "cat > \${RUNDIR}/varscan/vs_vep.indel.initial.input <<EOF\n";
    print VEP "varscan.vep.vcf = ./varscan.out.som_indel.gvip.vcf\n";
    print VEP "varscan.vep.output = ./varscan.out.som_indel.gvip.VEP.vcf\n";
    print VEP "varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VEP "varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VEP "varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VEP "varscan.vep.assembly = GRCh37\n";
    print VEP "EOF\n";
    print VEP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.snv.initial.input <<EOF\n";
    print VEP "strelka.vep.vcf = ./strelka.somatic.snv.strlk_pass.gvip.vcf\n";
    print VEP "strelka.vep.output = ./strelka.somatic.snv.strlk_pass.gvip.VEP.vcf\n";
    print VEP "strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VEP "strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VEP "strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VEP "strelka.vep.assembly = GRCh37\n";
    print VEP "EOF\n";
    print VEP "cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.indel.initial.input <<EOF\n";
    print VEP "strelka.vep.vcf = ./strelka.somatic.indel.strlk_pass.gvip.vcf\n";
    print VEP "strelka.vep.output = ./strelka.somatic.indel.strlk_pass.gvip.VEP.vcf\n";
    print VEP "strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print VEP "strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print VEP "strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print VEP "strelka.vep.assembly = GRCh37\n";
    print VEP "EOF\n";
	print VEP "cd \${RUNDIR}/varscan\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./vs_vep.snv.input >& ./vs_vep.snv.log\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./vs_vep.indel.input >& ./vs_vep.indel.log\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./vs_vep.snv.initial.input >& ./vs_vep.snv.initial.log\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./vs_vep.indel.initial.input >& ./vs_vep.indel.initial.log\n";
    print VEP "cd \${RUNDIR}/strelka/strelka_out/results\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./strelka_vep.snv.input >& ./strelka_vep.snv.log\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./strelka_vep.indel.input >& ./strelka_vep.indel.log\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./strelka_vep.snv.initial.input >& ./strelka_vep.snv.initial.log\n";
    print VEP "     ".$run_script_path."vep_annotator.pl ./strelka_vep.indel.initial.input >& ./strelka_vep.indel.initial.log\n";
	close VEP;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

}

sub bsub_parse_pindel {

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j8_parse_pindel".$sample_name.".sh";

    open(PP, ">$job_files_dir/$current_job_file") or die $!;
    print PP "#!/bin/bash\n";
    print PP "#BSUB -n 1\n";
    print PP "#BSUB -R \"rusage[mem=30000]\"","\n";
    print PP "#BSUB -M 30000000\n";
    print PP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print PP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print PP "#BSUB -J $current_job_file\n";
    print PP "#BSUB -q ding-lab\n";
    print PP "#BSUB -w \"$hold_job_file\"","\n";
    print PP "RUNDIR=".$sample_full_path."\n";
	print PP "cat > \${RUNDIR}/pindel/pindel_filter.input <<EOF\n";
	print PP "pindel.filter.pindel2vcf = $PINDEL_DIR/pindel2vcf\n";
	print PP "pindel.filter.variants_file = \${RUNDIR}/pindel/pindel.out.raw\n";
	print PP "pindel.filter.REF = $h37_REF\n";
	print PP "pindel.filter.date = 000000\n";
	print PP "pindel.filter.heterozyg_min_var_allele_freq = 0.2\n";
	print PP "pindel.filter.homozyg_min_var_allele_freq = 0.8\n";
	print PP "pindel.filter.mode = somatic\n";
	print PP "pindel.filter.apply_filter = true\n";
	print PP "pindel.filter.somatic.min_coverages = 10\n";
	print PP "pindel.filter.somatic.min_var_allele_freq = 0.10\n";
	print PP "pindel.filter.somatic.require_balanced_reads = \"true\"\n";
	print PP "pindel.filter.somatic.remove_complex_indels = \"true\"\n";
	print PP "pindel.filter.somatic.max_num_homopolymer_repeat_units = 6\n";
	print PP "EOF\n";
    print PP "cat > \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input <<EOF\n";
    print PP "pindel.dbsnp.indel.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar\n";
    print PP "pindel.dbsnp.indel.db = /gscmnt/gc3027/dinglab/medseq/cosmic/00-All.brief.pass.cosmic.vcf\n";
    print PP "pindel.dbsnp.indel.rawvcf = ./pindel.out.current_final.gvip.Somatic.vcf\n";
    print PP "pindel.dbsnp.indel.mode = filter\n";
    print PP "pindel.dbsnp.indel.passfile  = ./pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
    print PP "pindel.dbsnp.indel.dbsnpfile = ./pindel.out.current_final.gvip.dbsnp_present.vcf\n";
    print PP "EOF\n";
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
    print PP "     ".$run_script_path."dbsnp_filter.pl \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input\n";	
	print PP "cat > \${RUNDIR}/pindel/pindel_vep.input <<EOF\n";
	print PP "pindel.vep.vcf = ./pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
	print PP "pindel.vep.output = ./pindel.out.current_final.gvip.Somatic.VEP.vcf\n";
	print PP "pindel.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
	print PP "pindel.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
	print PP "pindel.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
	print PP "pindel.vep.assembly = GRCh37\n";
	print PP "EOF\n";
	print PP "     ".$run_script_path."vep_annotator.pl ./pindel_vep.input >& ./pindel_vep.log\n";  
	close PP;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
	}

sub bsub_merge_vcf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j9_merge_vcf.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    print MERGE "#BSUB -n 1\n";
    print MERGE "#BSUB -R \"rusage[mem=30000]\"","\n";
    print MERGE "#BSUB -M 30000000\n";
    print MERGE "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print MERGE "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print MERGE "#BSUB -J $current_job_file\n";
    print MERGE "#BSUB -q ding-lab\n";
    print MERGE "#BSUB -w \"$hold_job_file\"","\n";
    print MERGE "scr_t0=\`date \+\%s\`\n";
    print MERGE "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print MERGE "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print MERGE "myRUNDIR=".$sample_full_path."/varscan\n";
    print MERGE "STATUSDIR=".$sample_full_path."/status\n";
    print MERGE "RUNDIR=".$sample_full_path."\n";
    #print VEP "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print MERGE "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print MERGE "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print MERGE "export JAVA_OPTS=\"-Xmx2g\"\n";
    print MERGE "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
	print MERGE "STRELKA_VCF="."\${RUNDIR}/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf\n";
	print MERGE "VARSCAN_VCF="."\${RUNDIR}/varscan/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf\n";
	print MERGE "PINDEL_VCF="."\${RUNDIR}/pindel/pindel.out.current_final.gvip.Somatic.vcf\n";
	print MERGE "VARSCAN_INDEL="."\${RUNDIR}/varscan/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf\n";
	print MERGE "MERGER_OUT="."\${RUNDIR}/merged.vcf\n";
    print MERGE "cat > \${RUNDIR}/vep.merged.input <<EOF\n";
    print MERGE "merged.vep.vcf = ./merged.vcf\n"; 
    print MERGE "merged.vep.output = ./merged.VEP.vcf\n";
    print MERGE "merged.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print MERGE "merged.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print MERGE "merged.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print MERGE "merged.vep.assembly = GRCh37\n";
    print MERGE "EOF\n";
	print MERGE "java \${JAVA_OPTS} -jar $gatk -R $h37_REF -T CombineVariants -o \${MERGER_OUT} --variant:varscan \${VARSCAN_VCF} --variant:strelka \${STRELKA_VCF} --variant:varindel \${VARSCAN_INDEL} --variant:pindel \${PINDEL_VCF} -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,pindel,varindel\n"; 
	print MERGE "cd \${RUNDIR}\n";
	print MERGE "     ".$run_script_path."vep_annotator.pl ./vep.merged.input >&./vep.merged.log\n";	
	close MERGE;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
	}

sub bsub_vcf_2_maf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j10_vcf_2_maf.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    open(MAF, ">$job_files_dir/$current_job_file") or die $!;
    print MAF "#!/bin/bash\n";
    print MAF "#BSUB -n 1\n";
    print MAF "#BSUB -R \"rusage[mem=30000]\"","\n";
    print MAF "#BSUB -M 30000000\n";
    print MAF "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print MAF "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print MAF "#BSUB -J $current_job_file\n";
    print MAF "#BSUB -q ding-lab\n";
    print MAF "#BSUB -w \"$hold_job_file\"","\n";
    print MAF "F_VCF_1=".$sample_full_path."/merged.vcf\n";
	print MAF "F_VCF_2=".$sample_full_path."/".$sample_name.".vcf\n";
    print MAF "F_VEP_1=".$sample_full_path."/merged.VEP.vcf\n";
    print MAF "F_VEP_2=".$sample_full_path."/".$sample_name.".vep.vcf\n";
	print MAF "F_maf=".$sample_full_path."/".$sample_name.".maf\n";
	print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
	print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
	print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf	\${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $h37_REF --filter-vcf $f_exac\n";
    close MAF;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

	}
 
