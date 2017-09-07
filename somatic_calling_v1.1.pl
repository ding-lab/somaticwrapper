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
require('src/bsub_strelka.pl');
require("src/bsub_varscan.pl");
require("src/bsub_parse_strelka.pl");
require("src/bsub_parse_varscan.pl");

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
# everything else below should be automated
my $working_name= (split(/\//,$run_dir))[-2];
my $HOME1="/usr/local/somaticwrapper/runtime";
#store job files here
if (! -d $HOME1."/tmpsomatic2") {
    `mkdir -p $HOME1"/tmpsomatic2"`;
}
my $job_files_dir = $HOME1."/tmpsomatic2";
#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR_SOMATIC2") {
    `mkdir -p $HOME1"/LSF_DIR_SOMATIC2"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR_SOMATIC2";

# Distinguising between location of modules of somatic wrapper and GenomeVIP
my $script_dir="/usr/local/somaticwrapper/GenomeVIP";
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
#my $REF="/data/A_Reference/GRCh37-lite.fa";  # This does not work with strelka demo data because wrong reference
my $REF="/data/A_Reference/demo20.fa";  # Default strelka reference
my $f_exac="/gscmnt/gc2741/ding/qgao/tools/vcf2maf-1.6.11/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
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
                    bsub_strelka($sample_name, $sample_full_path, $job_files_dir, $bsub, $STRELKA_DIR, $REF);
                } elsif ($step_number == 2) {
                    bsub_varscan($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF);
                } 
#elsif ($step_number == 3) {
#                   &bsub_mutect(1);
#               }
                elsif ($step_number == 3) {
                    bsub_parse_strelka($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $script_dir);
                } 
                elsif ($step_number == 4) {
                    bsub_parse_varscan($sample_name, $sample_full_path, $job_files_dir, $bsub, $REF, $perl, $script_dir);
                } elsif ($step_number == 5) {
                    &bsub_pindel(1);
                }elsif ($step_number == 6) {
                    &bsub_vep(1);
                }elsif ($step_number == 7) {
                    &bsub_parse_pindel(1);
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








sub bsub_pindel{
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j5_pindel".$sample_name.".sh";  
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
    print PINDEL "$pindel -T 4 -f $REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"." -m 6 -w 1 -J $f_centromere\n";
    close PINDEL;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );   
}

sub bsub_vep{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j6_vep".$sample_name.".sh";
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
    print VEP "#BSUB -q long\n";
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
    print VEP "$perl $script_dir/vep_annotator.pl ./vs_vep.snv.input >& ./vs_vep.snv.log\n";
    print VEP "$perl $script_dir/vep_annotator.pl ./vs_vep.indel.input >& ./vs_vep.indel.log\n";
    print VEP "$perl $script_dir/vep_annotator.pl ./vs_vep.snv.initial.input >& ./vs_vep.snv.initial.log\n";
    print VEP "$perl $script_dir/vep_annotator.pl ./vs_vep.indel.initial.input >& ./vs_vep.indel.initial.log\n";
    print VEP "cd \${RUNDIR}/strelka/strelka_out/results\n";
    print VEP "$perl $script_dir/vep_annotator.pl ./strelka_vep.snv.input >& ./strelka_vep.snv.log\n";
    print VEP "$perl $script_dir/vep_annotator.pl ./strelka_vep.indel.input >& ./strelka_vep.indel.log\n";
    print VEP "$perl $script_dir/vep_annotator.pl ./strelka_vep.snv.initial.input >& ./strelka_vep.snv.initial.log\n";
    print VEP "$perl $script_dir/vep_annotator.pl ./strelka_vep.indel.initial.input >& ./strelka_vep.indel.initial.log\n";
    close VEP;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

}

sub bsub_parse_pindel {

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j7_parse_pindel".$sample_name.".sh";

    open(PP, ">$job_files_dir/$current_job_file") or die $!;
    print PP "#!/bin/bash\n";
    print PP "#BSUB -n 1\n";
    print PP "#BSUB -R \"rusage[mem=30000]\"","\n";
    print PP "#BSUB -M 30000000\n";
    print PP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print PP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print PP "#BSUB -J $current_job_file\n";
    print PP "#BSUB -q long\n";
    print PP "#BSUB -w \"$hold_job_file\"","\n";
    print PP "RUNDIR=".$sample_full_path."\n";
    print PP "cat > \${RUNDIR}/pindel/pindel_filter.input <<EOF\n";
    print PP "pindel.filter.pindel2vcf = $PINDEL_DIR/pindel2vcf\n";
    print PP "pindel.filter.variants_file = \${RUNDIR}/pindel/pindel.out.raw\n";
    print PP "pindel.filter.REF = $REF\n";
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
    print PP "$perl $script_dir/pindel_filter.v0.5.pl ./pindel_filter.input\n"; 
    print PP 'pre_current_final=$pin_var_file.CvgVafStrand_pass.Homopolymer_pass.vcf'."\n";
    print PP 'for mytmp in $pin_var_file.CvgVafStrand_pass.vcf  $pre_current_final  ${pre_current_final/%pass.vcf/fail.vcf} ; do'."\n";
    print PP '$perl $script_dir/genomevip_label.pl Pindel ./$mytmp ./${mytmp/%vcf/gvip.vcf}'."\n";
    print PP "done\n";
    print PP 'current_final=${pin_var_file/%raw/current_final.gvip.Somatic.vcf}'."\n";
    print PP 'cat ./${pre_current_final/%vcf/gvip.vcf} > ./$current_final'."\n";
    print PP "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print PP "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print PP "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print PP "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
    print PP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print PP "else\n";
    print PP "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print PP "fi\n";
    print PP "$perl $script_dir/dbsnp_filter.pl \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input\n";    
    print PP "cat > \${RUNDIR}/pindel/pindel_vep.input <<EOF\n";
    print PP "pindel.vep.vcf = ./pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
    print PP "pindel.vep.output = ./pindel.out.current_final.gvip.dbsnp_pass.VEP.vcf\n";
    print PP "pindel.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print PP "pindel.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print PP "pindel.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print PP "pindel.vep.assembly = GRCh37\n";
    print PP "EOF\n";
    print PP "$perl $script_dir/vep_annotator.pl ./pindel_vep.input >& ./pindel_vep.log\n";  
    close PP;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}

sub bsub_merge_vcf{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j8_merge_vcf.".$sample_name.".sh";
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
    print MERGE "#BSUB -q long\n";
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
    print MERGE "PINDEL_VCF="."\${RUNDIR}/pindel/pindel.out.current_final.gvip.dbsnp_pass.vcf\n";
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
    print MERGE "java \${JAVA_OPTS} -jar $gatk -R $REF -T CombineVariants -o \${MERGER_OUT} --variant:varscan \${VARSCAN_VCF} --variant:strelka \${STRELKA_VCF} --variant:varindel \${VARSCAN_INDEL} --variant:pindel \${PINDEL_VCF} -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,pindel,varindel\n"; 
    print MERGE "cd \${RUNDIR}\n";
    print MERGE "$perl $script_dir/vep_annotator.pl ./vep.merged.input >&./vep.merged.log\n";    
    close MERGE;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}

sub bsub_vcf_2_maf{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j9_vcf_2_maf.".$sample_name.".sh";
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
    print MAF "#BSUB -q long\n";
    print MAF "#BSUB -w \"$hold_job_file\"","\n";
    print MAF "F_VCF_1=".$sample_full_path."/merged.vcf\n";
    print MAF "F_VCF_2=".$sample_full_path."/".$sample_name.".vcf\n";
    print MAF "F_VEP_1=".$sample_full_path."/merged.VEP.vcf\n";
    print MAF "F_VEP_2=".$sample_full_path."/".$sample_name.".vep.vcf\n";
    print MAF "F_maf=".$sample_full_path."/".$sample_name.".maf\n";
    print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
    print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
    print MAF "$perl $script_dir/vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $REF --filter-vcf $f_exac\n";
    close MAF;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

}

