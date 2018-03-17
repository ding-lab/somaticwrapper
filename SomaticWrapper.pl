######### Song Cao and Matt Wyczalkowski ###########
## pipeline for somatic variant calling ##

# TODO: implement full command line argument parsing instead of configuration file
# see: https://perlmaven.com/how-to-process-command-line-arguments-in-perl

# TODO: how are paths handled in CGC?  Do we need so pass output paths individually?

# sample_name -> run_name

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

#use POSIX;
my $version = 3.0;

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
This script will evaluate variants for WGS and WXS data
Pipeline version: $version
Usage: perl $0 [options] run_name step_number 

run_name is unique identifier of this analysis (previously called sample_name)
step_number executes given step of pipeline:
* [1 or run_strelka]  Run streka
* [2 or run_varscan]  Run Varscan
* [3 or parse_strelka]  Parse streka result
* [4 or parse_varscan]  Parse VarScan result
* [5 or run_pindel]  Run Pindel
* [7 or parse_pindel]  Parse Pindel
* [8 or merge_vcf]  Merge vcf files using local VEP cache
* [10 or run_vep]  Run VEP annotation

Optional configuration file parameters [defaults]
    --tumor_bam s:  path to tumor BAM.  Required for all runs
    --normal_bam s: path to normal BAM.  Required for all runs
    --reference_fasta s: path to reference
    --assembly s: either "GRCh37" or "GRCh38", used for VEP [GRCh37]
    --reference_dict s: path to reference dict file.  Default is reference_fasta with ".dict" appended
    --sw_dir s: Somatic Wrapper installation directory [/usr/local/somaticwrapper]
    --sw_data s: Somatic Wrapper analysis results directory [/data/data]
            Per-sample analysis directory is sw_data/run_name
    --use_vep_db : if defined, use online VEP database lookups ("db mode") [false]
          db mode a) uses online database (so cache isn't installed) b) does not use tmp files
          It is meant to be used for testing and lightweight applications.  Use the cache (default)
          for better performance.
          See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
    --vep_cache_dir s: VEP cache directory, if not doing online VEP db lookups.  [/data/D_VEP]
    --output_vep : if defined, write final annotated merged file in VEP rather than VCF format [false]
    --annotate_intermediate : if defined, annotate intermediate files with VEP [false]
    --strelka_config s: path to strelka.ini file, required for strelka run
    --varscan_config s: path to varscan.ini file, required for varscan run
    --pindel_config s: path to pindel.ini file, required for pindel parsing
    --centromere_bed s: path to BED file describing centromere regions to exclude for pindel analysis.  See C_Centromeres for discussion
    --gatk_jar s: path to GATK Jar file.  [/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar]
    --perl s: path to PERL executable.  [/usr/bin/perl]
    --strelka_dir s: path to strelka installation dir.  [/usr/local/strelka]
    --vep_cmd s: path to ensembl vep executable.  [/usr/local/ensembl-vep/vep]
    --pindel_dir s: path to Pindel installation dir.  [/usr/local/pindel]
    --snpsift_jar s: [/usr/local/snpEff/SnpSift.jar]
    --varscan_jar s: [usr/local/VarScan.jar]
    --dbsnp_db s: database for dbSNP filtering [none]
OUT

# OLD:
# my $centromere_bed="$sw_dir/image.setup/C_Centromeres/pindel-centromere-exclude.bed";

# Argument parsing reference: http://perldoc.perl.org/Getopt/Long.html
# https://perlmaven.com/how-to-process-command-line-arguments-in-perl
my $tumor_bam;
my $normal_bam;
my $assembly="GRCh37";
my $reference_fasta;
my $reference_dict;  # default mapping occurs after reference_fasta known
my $sw_dir = "/usr/local/somaticwrapper";
my $sw_data = "/data/data";
my $use_vep_db; 
my $vep_cache_dir = "/data/D_VEP";
my $output_vep;
my $annotate_intermediate;
my $strelka_config; 
my $varscan_config; 
my $pindel_config; 
my $centromere_bed; 
my $gatk_jar = "/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar";
my $perl = "/usr/bin/perl";
my $strelka_dir = "/usr/local/strelka";
my $vep_cmd = "/usr/local/ensembl-vep/vep";
my $pindel_dir = "/usr/local/pindel";
my $snpsift_jar = "/usr/local/snpEff/SnpSift.jar";
my $varscan_jar = "usr/local/VarScan.jar";
my $dbsnp_db = "none";

GetOptions(
    'tumor_bam=s' => \$tumor_bam,
    'normal_bam=s' => \$normal_bam,
    'reference_fasta=s' => \$reference_fasta,
    'assembly=s' => \$assembly,
    'reference_dict=s' => \$reference_dict,
    'sw_dir=s' => \$sw_dir,
    'sw_data=s' => \$sw_data,
    'use_vep_db' => \$use_vep_db,
    'vep_cache_dir=s' => \$vep_cache_dir,
    'output_vep' => \$output_vep,
    'annotate_intermediate' => \$annotate_intermediate,
    'strelka_config=s' => \$strelka_config,
    'varscan_config=s' => \$varscan_config,
    'pindel_config=s' => \$pindel_config,
    'centromere_bed=s' => \$centromere_bed,
    'gatk_jar=s' => \$gatk_jar,
    'perl=s' => \$perl,
    'strelka_dir=s' => \$strelka_dir,
    'vep_cmd=s' => \$vep_cmd,
    'pindel_dir=s' => \$pindel_dir,
    'snpsift_jar=s' => \$snpsift_jar,
    'varscan_jar=s' => \$varscan_jar
    'dbsnp_db=s' => \$dbsnp_db
) or die "Error parsing command line args.\n$usage\n";

die $usage unless @ARGV >= 2;
my ( $run_name, $step_number ) = @ARGV;

if ( not $reference_dict ) $reference_dict = "$reference_fasta.dict";




die("tumor_bam undefined in $config_file\n") unless $tumor_bam;
die("normal_bam undefined in $config_file\n") unless $normal_bam;
die("assembly undefined in $config_file\n") unless $assembly;
die("reference_fasta undefined in $config_file\n") unless $reference_fasta;
die("reference_dict undefined in $config_file\n") unless $reference_dict;
die("run_name undefined in $config_file\n") unless $run_name;

# Distinguising between location of modules of somatic wrapper and GenomeVIP
# GenomeVIP is not distributed separately so hard code the path
my $gvip_dir="$sw_dir/GenomeVIP";

#begin to process each sample
my $sample_full_path = $sw_data."/".$run_name;

# automatically generated scripts in runtime
my $job_files_dir="$sample_full_path/runtime";
system("mkdir -p $job_files_dir");

print("Using reference $ref\n");
print("SomaticWrapper dir: $sw_dir \n");
print("Analysis dir: $sample_full_path\n");
print("Run script dir: $job_files_dir\n");

print "\nSubmitting jobs for the sample ",$run_name, "..." "\n";

if (($step_number eq '1') || ($step_number eq 'run_strelka')) {
    die("strelka_config undefined in $config_file\n") unless exists $paras{'strelka_config'};
#    my $strelka_config = "/usr/local/somaticwrapper/config/strelka.ini";
    run_strelka($tumor_bam, $normal_bam, $run_name, $sample_full_path, $job_files_dir, $bsub, $strelka_dir, $ref, $paras{'strelka_config'});
} elsif (($step_number eq '2') || ($step_number eq 'run_varscan')) {
    die("varscan_config undefined in $config_file\n") unless exists $paras{'varscan_config'};
    run_varscan($tumor_bam, $normal_bam, $run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $paras{'varscan_config'});
} elsif (($step_number eq '3') || ($step_number eq 'parse_strelka')) {
    parse_strelka($run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $ref_dict, $perl, $gvip_dir, $dbsnp_db, $snpsift_jar);
} elsif (($step_number eq '4') || ($step_number eq 'parse_varscan')) {
    parse_varscan($run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir, $dbsnp_db, $snpsift_jar, $varscan_jar);
} elsif (($step_number eq '5') || ($step_number eq 'run_pindel')) {
    run_pindel($tumor_bam, $normal_bam, $run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $pindel_dir, $centromere_bed);
} elsif (($step_number eq '7') || ($step_number eq 'parse_pindel')) {
    die("pindel_config undefined in $config_file\n") unless exists $paras{'pindel_config'};
    parse_pindel($run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir, $vep_cmd, $pindel_dir, $dbsnp_db, $snpsift_jar, $paras{'pindel_config'});
} elsif (($step_number eq '8') || ($step_number eq 'merge_vcf')) {
    merge_vcf($run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir, $vep_cmd, $gatk_jar, $use_vep_db, $output_vep, $assembly, $vep_cache_dir);
} elsif (($step_number eq '9') || ($step_number eq 'vcf2maf')) {
    die("vcf_2_maf() disabled while ExAC CRCh38 issues resolved.\n");
    vcf_2_maf($run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir);
} elsif (($step_number eq '10') || ($step_number eq 'run_vep')) {
    print("annotate_intermediate = $annotate_intermediate\n");
    run_vep($run_name, $sample_full_path, $job_files_dir, $bsub, $ref, $gvip_dir, $vep_cmd, $assembly, $vep_cache_dir, $use_vep_db, $output_vep, $annotate_intermediate);
} else {
    die("Unknown step number $step_number\n");
}
