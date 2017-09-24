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
my $version = 3.0;
#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal_color = "\e[0m";
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
This script will process evaluate variants for WGS and WXS data
Pipeline version: $version
$yellow     Usage: perl $0 run_folder step_number config_file [config_file_2] $normal_color

run_dir = full path of the folder holding analysis resuls
            Note, per-sample analysis directory is run_dir/sample_name
step_number run this pipeline step by step. (running the whole pipeline if step number is 0)
config_file Input configuration file.  See below for format
config_file Optional secondary configuration file, any parameters here override configuration previous configuration

$green       [1]  Run streka
$red         [2]  Run Varscan
$yellow      [3]  Parse streka result
$purple      [4]  Parse VarScan result
$cyan        [5]  Run Pindel
$gray        [6]  Run VEP annotation
$gray        [7]  Parse Pindel
$gray        [8]  Merge vcf files using local VEP cache
$gray        [8b] Merge vcf files using VEP db queries 
$gray        [9] generate maf file 
$normal_color

Input File configuration

Format:
    key = value

Required configuration file keys
    tumor_bam
    normal_bam
    reference_fasta
    reference_dict
    sample_name

Optional configuration file parameters
    sw_dir - Somatic Wrapper installation directory
    usedb - whether to use online VEP database lookups (1 for true)
    submit_cmd - command to initiate execution of generated script.  Default value 'bash', can set as 'cat' to allow step-by-step execution for debugging

OUT

die $usage unless @ARGV == 3;
my ( $run_dir, $step_number, $config_file, $config_file2 ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {  # ?
    $run_dir = $1;
}

print("Reading configuration file $config_file\n");

# get paras from config file
# for a "key = value" pair of "xxx.yyy.zzz = foo", generates entry $params{'zzz'}='foo'
open(CONFIG, $config_file);
my (%paras);
map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<CONFIG>);
close(CONFIG);

# Goal of an optional secondary configuration file is to allow for local configuration changes which are not tracked by git.
# An example is to have a different sw_dir on MGI, but don't want MGI-specific changes propagated to main branch
# check if parameter is defined, and if file exists
if (defined $config_file2) {
    if (-e $config_file2) {
        print("Reading secondary configuration file $config_file2\n");
        open(CONFIG2, $config_file2);
        map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<CONFIG2>);
        close(CONFIG2);
    } else {
        print("Optional configuration file $config_file2 does not exist.  Continuing. \n");
    }
}

map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;



# Data and software configuration

# Data

# Required configuration file parameters
# tumor_bam
# normal_bam
# reference_fasta
# reference_dict
# sample_name

# Optional configuration file parameters
# sw_dir - Somatic Wrapper installation directory.  Default is /usr/local/somaticwrapper, modified for MGI installation
# usedb - whether to use online VEP database lookups.  If value is 0 or undefined, default to cache (which requires installation)
    # more detail from GenomeVIP/vep_annotator.pl:
    # db mode 1) uses online database (so cache isn't installed) 2) does not use tmp files
    # It is meant to be used for testing and lightweight applications.  Use the cache for
    # better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 


# default convention:
#$sample_full_path = $run_dir."/".$sample_name;
#my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
#my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
die("tumor_bam undefined in $config_file\n") unless exists $paras{'tumor_bam'};
my $tumor_bam = $paras{'tumor_bam'};

die("normal_bam undefined in $config_file\n") unless exists $paras{'normal_bam'};
my $normal_bam = $paras{'normal_bam'};

die("reference_fasta undefined in $config_file\n") unless exists $paras{'reference_fasta'};
my $REF = $paras{'reference_fasta'};

die("reference_dict undefined in $config_file\n") unless exists $paras{'reference_dict'};
my $ref_dict = $paras{'reference_dict'};

die("sample_name undefined in $config_file\n") unless exists $paras{'sample_name'};
my $sample_name = $paras{'sample_name'};

# This is the default somatic wrapper installation directory
my $sw_dir="/usr/local/somaticwrapper";
if (exists $paras{'sw_dir'} ) {
    $sw_dir=$paras{'sw_dir'};
}

# This is the default somatic wrapper installation directory
my $usedb=0;
if (exists $paras{'usedb'} ) {
    $usedb=$paras{'usedb'};
}

# parameter submit_cmd will typically be "bash" to execute entire script.  Set it to "cat" for debugging, "bsub" to submit to LSF 
my $bsub = "bash";
if (exists $paras{'submit_cmd'} ) {
    $bsub=$paras{'submit_cmd'};
}



# TODO: allow custom dbSnP/COSMIC filter, so that Strelka Demo filter does not clobber real filter
# this will be an optional parameter with default value of dbsnp.noCOSMIC.vcf.gz

# Define where centromere definion file is for pindel processing.  See C_Centromeres for discussion
# This should really live in the data directory
# This should be defined in configuration file for pindel step
my $f_centromere="$sw_dir/image.setup/C_Centromeres/pindel-centromere-exclude.bed";

# Distinguising between location of modules of somatic wrapper and GenomeVIP
my $gvip_dir="$sw_dir/GenomeVIP";
my $perl = "/usr/bin/perl";
#my $STRELKA_DIR="/usr/local/strelka-2.7.1.centos5_x86_64/bin";
# using older version of strelka
my $STRELKA_DIR="/usr/local/strelka";
my $vep_cmd="/usr/local/ensembl-vep/vep";
my $pindel_dir="/usr/local/pindel";
my $gatk="/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar";

#begin to process each sample
my $sample_full_path = $run_dir."/".$sample_name;

# automatically generated scripts in runtime
my $job_files_dir="$sample_full_path/runtime";
system("mkdir -p $job_files_dir");

print("Using reference $REF\n");
print("SomaticWrapper dir: $sw_dir \n");
print("Analysis dir: $sample_full_path\n");
print("Run script dir: $job_files_dir\n");

if (-d $sample_full_path) { # is a full path directory containing sample analysis
    print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal_color, "\n";

    if ($step_number eq '1') {
        run_strelka($tumor_bam, $normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $STRELKA_DIR, $REF);
    } elsif ($step_number eq '2') {
        run_varscan($tumor_bam, $normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $REF);
    } elsif ($step_number eq '3') {
        parse_strelka($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $ref_dict, $perl, $gvip_dir);
    } elsif ($step_number eq '4') {
        parse_varscan($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $gvip_dir);
    } elsif ($step_number eq '5') {
        run_pindel($tumor_bam, $normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $pindel_dir, $f_centromere);
    } elsif ($step_number eq '6') {
        warn("run_vep() is ignored in pipeline, so output of this step is discarded.  Continuing anyway.\n\n");
        run_vep($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $gvip_dir, $vep_cmd);
    } elsif ($step_number eq '7') {
        parse_pindel($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $gvip_dir, $vep_cmd, $pindel_dir);
    } elsif ($step_number eq '8') {
        merge_vcf($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $gvip_dir, $vep_cmd, $gatk, $usedb);
    } elsif ($step_number eq '9') {
        die("vcf_2_maf() disabled while ExAC CRCh38 issues resolved.\n");
        vcf_2_maf($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $gvip_dir);
    } else {
        die("Unknown step number $step_number\n");
    }
}

