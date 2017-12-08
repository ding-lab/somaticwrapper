######### Song Cao and Matt Wyczalkowski ###########
## pipeline for somatic variant calling ##

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

run_dir = full path of the folder holding analysis results
            Note, per-sample analysis directory is run_dir/sample_name
step_number run this pipeline step by step. (running the whole pipeline if step number is 0)
config_file Input configuration file.  See below for format
config_file_2 Optional secondary configuration file, any parameters here override configuration previous configuration

$green       [1 or run_strelka]  Run streka
$red         [2 or run_varscan]  Run Varscan
$yellow      [3 or parse_strelka]  Parse streka result
$purple      [4 or parse_varscan]  Parse VarScan result
$cyan        [5 or run_pindel]  Run Pindel
$gray        [7 or parse_pindel]  Parse Pindel
$gray        [8 or merge_vcf]  Merge vcf files using local VEP cache
$gray        [9 or vcf2maf] generate maf file 
$gray        [10 or run_vep]  Run VEP annotation
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
    assembly - GRCh37 or GRCh38

Optional configuration file parameters
    sw_dir - Somatic Wrapper installation directory.  Default is /usr/local/somaticwrapper
    use_vep_db - whether to use online VEP database lookups (1 for true)
          db mode a) uses online database (so cache isn't installed) b) does not use tmp files
          It is meant to be used for testing and lightweight applications.  Use the cache for
          better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
    vep_cache_dir - VEP cache directory, if not doing online VEP db lookups.  Default is "/data/D_VEP"
    submit_cmd - command to initiate execution of generated script.  Default value 'bash', can set as 'cat' to allow step-by-step execution for debugging
    output_vep - write final annotated merged file in VEP rather than VCF format
    annotate_intermediate - VEP-annotate intermediate output files
    strelka_config - path to strelka.ini file, required for strelka run
    varscan_config - path to varscan.ini file, required for varscan run
    centromere_bed - path to BED file describing centromere regions to exclude for pindel analysis.  
        Default: $sw_dir/image.setup/C_Centromeres/pindel-centromere-exclude.bed
    gatk - path to GATK Jar file.  Default: /usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
    perl - path to PERL executable.  Default: /usr/bin/perl
    strelka_dir - path to strelka installation dir.  Default: /usr/local/strelka
    vep_cmd - path to ensembl vep executable.  Default: /usr/local/ensembl-vep/vep
    pindel_dir - path to Pindel installation dir.  Default: /usr/local/pindel

OUT

die $usage unless @ARGV == 3;
my ( $run_dir, $step_number, $config_file, $config_file2 ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {  # ?
    $run_dir = $1;
}

print("Reading configuration file $config_file\n");

# get paras from config file
# for a "key = value" pair of "xxx.yyy.zzz = foo", generates entry $params{'zzz'}='foo'
open(CONFIG, $config_file) or die "Could not open file '$config_file' $!";
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

die("tumor_bam undefined in $config_file\n") unless exists $paras{'tumor_bam'};
my $tumor_bam = $paras{'tumor_bam'};

die("normal_bam undefined in $config_file\n") unless exists $paras{'normal_bam'};
my $normal_bam = $paras{'normal_bam'};

die("assembly undefined in $config_file\n") unless exists $paras{'assembly'};
my $assembly = $paras{'assembly'};

die("reference_fasta undefined in $config_file\n") unless exists $paras{'reference_fasta'};
my $ref = $paras{'reference_fasta'};

die("reference_dict undefined in $config_file\n") unless exists $paras{'reference_dict'};
my $ref_dict = $paras{'reference_dict'};

die("sample_name undefined in $config_file\n") unless exists $paras{'sample_name'};
my $sample_name = $paras{'sample_name'};

# This is the default somatic wrapper installation directory
my $sw_dir="/usr/local/somaticwrapper";
if (exists $paras{'sw_dir'} ) {
    $sw_dir=$paras{'sw_dir'};
}

my $use_vep_db=0;
if (exists $paras{'use_vep_db'} ) {
    $use_vep_db=$paras{'use_vep_db'};
}

# parameter submit_cmd will typically be "bash" to execute entire script.  Set it to "cat" for debugging, "bsub" to submit to LSF 
my $bsub = "bash";
if (exists $paras{'submit_cmd'} ) {
    $bsub=$paras{'submit_cmd'};
}

my $dbsnp_db = "none";  # TODO: if dbsnp_db not defined, either give an error or skip dbSNP filtering
if (exists $paras{'dbsnp_db'} ) {
    $dbsnp_db=$paras{'dbsnp_db'};
}

my $vep_cache_dir = "/data/D_VEP";
if (exists $paras{'vep_cache_dir'} ) {
    $vep_cache_dir=$paras{'vep_cache_dir'};
}

my $output_vep = 0;
if (exists $paras{'output_vep'} ) {
    $output_vep=$paras{'output_vep'};
}

my $annotate_intermediate=0;
if (exists $paras{'annotate_intermediate'} ) {
    $annotate_intermediate=$paras{'annotate_intermediate'};
}

# Add params below to optional params with defaults

# Define where centromere definion file is for pindel processing.  See C_Centromeres for discussion
# This should really live in the data directory
my $centromere_bed="$sw_dir/image.setup/C_Centromeres/pindel-centromere-exclude.bed";
if (exists $paras{'centromere_bed'} ) {
    $centromere_bed=$paras{'centromere_bed'};
}

my $gatk="/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar";
if (exists $paras{'gatk'} ) {
    $gatk=$paras{'gatk'};
}

my $perl = "/usr/bin/perl";
if (exists $paras{'perl'} ) {
    $perl=$paras{'perl'};
}

#my $strelka_dir="/usr/local/strelka-2.7.1.centos5_x86_64/bin";
# using older version of strelka
my $strelka_dir="/usr/local/strelka";
if (exists $paras{'strelka_dir'} ) {
    $strelka_dir=$paras{'strelka_dir'};
}

my $vep_cmd="/usr/local/ensembl-vep/vep";
if (exists $paras{'vep_cmd'} ) {
    $vep_cmd=$paras{'vep_cmd'};
}

my $pindel_dir="/usr/local/pindel";
if (exists $paras{'pindel_dir'} ) {
    $pindel_dir=$paras{'pindel_dir'};
}

# Distinguising between location of modules of somatic wrapper and GenomeVIP
# GenomeVIP is not distributed separately so hard code the path
my $gvip_dir="$sw_dir/GenomeVIP";

#begin to process each sample
my $sample_full_path = $run_dir."/".$sample_name;

# automatically generated scripts in runtime
my $job_files_dir="$sample_full_path/runtime";
system("mkdir -p $job_files_dir");

print("Using reference $ref\n");
print("SomaticWrapper dir: $sw_dir \n");
print("Analysis dir: $sample_full_path\n");
print("Run script dir: $job_files_dir\n");

if (-d $sample_full_path) { # is a full path directory containing sample analysis
    print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal_color, "\n";

    if (($step_number eq '1') || ($step_number eq 'run_strelka')) {
        die("strelka_config undefined in $config_file\n") unless exists $paras{'strelka_config'};
    #    my $strelka_config = "/usr/local/somaticwrapper/config/strelka.ini";
        run_strelka($tumor_bam, $normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $strelka_dir, $ref, $paras{'strelka_config'});
    } elsif (($step_number eq '2') || ($step_number eq 'run_varscan')) {
        die("varscan_config undefined in $config_file\n") unless exists $paras{'varscan_config'};
        run_varscan($tumor_bam, $normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $paras{'varscan_config'});
    } elsif (($step_number eq '3') || ($step_number eq 'parse_strelka')) {
        parse_strelka($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $ref_dict, $perl, $gvip_dir, $dbsnp_db);
    } elsif (($step_number eq '4') || ($step_number eq 'parse_varscan')) {
        parse_varscan($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir, $dbsnp_db);
    } elsif (($step_number eq '5') || ($step_number eq 'run_pindel')) {
        run_pindel($tumor_bam, $normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $pindel_dir, $centromere_bed);
    } elsif (($step_number eq '7') || ($step_number eq 'parse_pindel')) {
        parse_pindel($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir, $vep_cmd, $pindel_dir, $dbsnp_db);
    } elsif (($step_number eq '8') || ($step_number eq 'merge_vcf')) {
        merge_vcf($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir, $vep_cmd, $gatk, $use_vep_db, $output_vep, $assembly, $vep_cache_dir);
    } elsif (($step_number eq '9') || ($step_number eq 'vcf2maf')) {
        die("vcf_2_maf() disabled while ExAC CRCh38 issues resolved.\n");
        vcf_2_maf($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir);
    } elsif (($step_number eq '10') || ($step_number eq 'run_vep')) {
        print("annotate_intermediate = $annotate_intermediate\n");
        run_vep($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $gvip_dir, $vep_cmd, $assembly, $vep_cache_dir, $use_vep_db, $output_vep, $annotate_intermediate);
    } else {
        die("Unknown step number $step_number\n");
    }
}

