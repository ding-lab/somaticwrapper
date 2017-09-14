######### Song Cao###########
## pipeline for somatic variant callings ##
#   somatic_variant_callings.pl #
### updated date: 04/05/2017 ###
### updated date: 04/18/2017 ###
### add vcf2maf.pl ###
### 07/14/2017 ##

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

# submodule information
# https://stackoverflow.com/questions/1712016/how-do-i-include-functions-from-another-file-in-my-perl-script
require('src/run_strelka.pl');
require("src/run_varscan.pl");
require("src/parse_strelka.pl");
require("src/parse_varscan.pl");
require("src/run_pindel.pl");
require("src/parse_pindel.pl");
require("src/run_vep.pl");
require("src/merge_vcf.pl");
require("src/vcf_2_maf.pl");

(my $usage = <<OUT) =~ s/\t+//g;
This script will process rna-seq data for TCGA samples. 
Pipeline version: $version
$yellow     Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$green       [1]  Run streka
$red         [2]  Run Varscan
$yellow      [3]  Parse streka result
$purple      [4]  Parse VarScan result
$cyan        [5]  Run Pindel
$gray        [6]  Run VEP annotation
$gray        [7]  Parse Pindel
$gray        [8]  Merge vcf files  
$gray        [9] generate maf file 
$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}
die $usage unless ($step_number >=0)&&(($step_number <= 10));

my $working_name= (split(/\//,$run_dir))[-2];

my $sw_dir="/usr/local/somaticwrapper";

# Distinguising between location of modules of somatic wrapper and GenomeVIP
my $gvip_dir="$sw_dir/GenomeVIP";

# automatically generated scripts below
my $job_files_dir="$sw_dir/runtime";
system("mkdir -p $job_files_dir");

# Define where centromere definion file is.  See C_Centromeres for discussion
my $f_centromere="$sw_dir/C_Centromeres/pindel-centromere-exclude.bed";

my $perl = "/usr/bin/perl";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub = "bash"; # or bsub 
my $sample_full_path = "";
my $sample_name = "";

#my $STRELKA_DIR="/usr/local/strelka-2.7.1.centos5_x86_64/bin";
# using older version of strelka
my $STRELKA_DIR="/usr/local/strelka";


# Note that VEP that had been used was v81, which used the script variant_effect_predictor.pl
# Newer versions of VEP use vep.pl as the command.  It remains to be seen what differences there are
# old: VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
my $vep_cmd="/usr/local/ensembl-vep/vep";

# For PINDEL testing, use the complete reference.
#my $REF="/data/A_Reference/GRCh37-lite.fa";  # This does not work with strelka demo data because wrong reference
my $REF="/data/A_Reference/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna";   # This is GRCh37 with 'chrN' naming
#my $REF="/data/A_Reference/demo20.fa";  # Strelka reference.  To use with vep, it has to have standard chrom names
print("Using reference $REF\n");

my $f_exac="/gscmnt/gc2741/ding/qgao/tools/vcf2maf-1.6.11/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
my $pindel_dir="/usr/local/pindel";
my $gatk="/gscuser/scao/tools/GenomeAnalysisTK.jar";

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
#  UNUSED           bsub_strelka();
#                    &bsub_varscan();
#&bsub_mutect();
#                   &bsub_parse_strelka();
#                   &bsub_parse_varscan();
#                   &bsub_pindel();
#                   &bsub_vep();
#                   &bsub_parse_pindel();
#                   &bsub_merge_vcf();
#                   &bsub_vcf_2_maf();
# &bsub_pindel();   
                }
                elsif ($step_number == 1) {
                    run_strelka($sample_name, $sample_full_path, $job_files_dir, $bsub, $STRELKA_DIR, $REF);
                } elsif ($step_number == 2) {
                    run_varscan($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF);
                } 
#elsif ($step_number == 3) {
#                   &bsub_mutect(1);
#               }
                elsif ($step_number == 3) {
                    parse_strelka($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $gvip_dir);
                } 
                elsif ($step_number == 4) {
                    parse_varscan($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $gvip_dir);
                } elsif ($step_number == 5) {
                    run_pindel($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $pindel_dir, $sw_dir, $f_centromere);
                }elsif ($step_number == 6) {
                    run_vep($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $gvip_dir, $vep_cmd);
                }elsif ($step_number == 7) {
                    parse_pindel($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $gvip_dir, $vep_cmd, $pindel_dir);
                }elsif ($step_number == 8) {
                    &bsub_merge_vcf(1);
                }elsif ($step_number == 9) {
                    &bsub_vcf_2_maf(1);
                }

# how to run pindel #
# /gscmnt/gc8001/info/model_data/8ad63791a648435b95a4be23ac3a18a5/buildb8e7ade8c3244bd483c274edeceeb8b7/alignments/cfa9c81317e346f5aa458e5e9220f41e.bam 500 TCGA-EX-A1H5-01A-31D-A34H-09 #
#/gscmnt/gc8001/info/model_data/c777ae8a6fdb49b7bf400056d15a6788/build8057e74cc9f64a50a13c33763af10116/alignments/b127fcaaecc54abb833316fa3153e492.bam 500 TCGA-EX-A1H5-10A-01D-A200-09#    
# bsub -q long -M 16000000 -n 4 -R 'span[hosts=1] rusage[mem=16000]' -oo /gscmnt/gc2532/dinglab/scao/pindel/CESC/pindel-logs/TCGA-C5-A0TN.chr9.log /gscuser/kye/gc2532/projects/PCGP_rerun/pindel -T 4 -f /gscmnt/gc2532/dinglab/projects/PCGP_rerun/GRCh37-lite.fa -i /gscmnt/gc2532/dinglab/scao/pindel/CESC/pindel-configs/TCGA-C5-A0TN.config -o /gscmnt/gc2532/dinglab/scao/pindel/CESC/pindel-outputs/TCGA-C5-A0TN/TCGA-C5-A0TN.chr9 -c 9 -m 6 -w 1 -J /gscmnt/gc3015/dinglab/medseq/Jiayin_Germline_Project/PCGP/data/pindel-centromere-exclude.bed #

            }
        }
    }
}

#######################################################################
if ($step_number == 0) {
    print $green, "All jobs are submitted!\n",$normal;
}

exit;











