######### SomaticWrapper: tumor-only pipeline (compute1) ###########
##### contact: scao@wustl.edu ####

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# ---- dependency globals ----
our @J1_JIDS      = ();   # per-sample, many jobs (one per chr)
our $J2_JID       = '';   # per-sample, single job (filter)
our $J3_JID       = '';   # per-sample, single job (parse)
our @ALL_J4_JIDS  = ();   # collect all samples' j4 job IDs (for run-level deps)
our $J5_JID       = '';   # run-level report job
# --------------------------------

my $version = 3.0;

# color for logs / usage
my $yellow = "\e[33m";
my $green  = "\e[32m";
my $cyan   = "\e[36m";
my $normal = "\e[0m";

(my $usage = <<OUT) =~ s/\t+//g;
Somatic variant calling pipeline 
Pipeline version: $version

$yellow     Usage: perl \$0 --chr --step --sre --rdir --ref --log --q --exonic --groupname --users $normal

<rdir> = full path of the folder holding files for this sequence run (required)
<log>  = full path of the folder for saving log files (required)
<chr>  = whether reference has 'chr' prefix (1 yes / 0 no; default 1)
<sre>  = rerun: 1 yes / 0 no (default 0)
<step> = run a stage (required): 
$green [0] Submit all steps with dependencies
$green [1] Run Mutect2
$green [2] Filter Mutect2 result
$green [3] Parse Mutect2 result
$cyan  [4] Generate per-sample MAF
$cyan  [5] Generate merged run-level report$normal

Options:
<ref>       human reference fasta (required)
<q>         LSF queue; e.g. long (default), ding-lab, research-hpc
<groupname> job group name used as /<users>/<groupname> (required)
<users>     compute1 username used in job group path (required)
<exonic>    exonic output flag for reporting: 1 yes / 0 no (default 1)
OUT

# ---- CLI args ----
my $step_number   = -1;
my $status_rerun  = 0;
my $status_exonic = 1;
my $q_name        = 'long';
my $run_dir       = '';
my $log_dir       = '';
my $h38_REF       = '';
my $chr_status    = 1;
my $compute_username = '';
my $group_name       = '';
my $mincov_t      = 14;      # used by parser
my $minvaf        = 0.05;    # used by parser

GetOptions(
  "step=i"      => \$step_number,
  "chr=i"       => \$chr_status,
  "sre=i"       => \$status_rerun,
  "groupname=s" => \$group_name,
  "users=s"     => \$compute_username,		
  "exonic=i"    => \$status_exonic,
  "rdir=s"      => \$run_dir,
  "ref=s"       => \$h38_REF,
  "log=s"       => \$log_dir,
  "q=s"         => \$q_name,
) or die $usage;

if (
   $run_dir eq "" || $log_dir eq "" || $h38_REF eq "" ||
   $step_number < 0 || $group_name eq "" || $compute_username eq ""
) {
  print "wrong option\n$usage";
  exit 1;
}

print "run dir=$run_dir\n";
print "log dir=$log_dir\n";
print "step num=$step_number\n";
print "status rerun=$status_rerun\n";
print "status exonic=$status_exonic\n";
print "queue name=$q_name\n";

$run_dir =~ s/(.+)\/$/$1/;
die $usage unless ($step_number >= 0 && $step_number <= 5);

# ---- paths & dirs ----
my $HOME1 = $log_dir;

if (! -d $HOME1)                  { `mkdir -p $HOME1`; }
if (! -d "$HOME1/tmpsomatic")     { `mkdir -p $HOME1/tmpsomatic`; }
if (! -d "$HOME1/LSF_DIR_SOMATIC"){ `mkdir -p $HOME1/LSF_DIR_SOMATIC`; }

my $job_files_dir = "$HOME1/tmpsomatic";
my $lsf_file_dir  = "$HOME1/LSF_DIR_SOMATIC";

# script path prefix for helper scripts (perl)
my $pwd = `echo \$PWD`; chomp $pwd;
my $run_script_path = "/usr/bin/perl $pwd/";

# java & tools actually used below
my $java_dir  = "/storage1/fs1/songcao/Active/Software/jre1.8.0_121";
my $java_bin  = "/storage1/fs1/songcao/Active/Software/jre1.8.0_121/bin/java";
my $GATK      = "/storage1/fs1/dinglab/Active/Projects/austins2/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar";
my $picardexe = "/storage1/fs1/songcao/Active/Software/picard/picard.jar";
my $snpsift   = "/storage1/fs1/songcao/Active/Software/snpEff/20150522/SnpSift.jar";
my $vepannot  = "/storage1/fs1/dinglab/Active/Projects/scao/gitshared/ensembl-vep/vep";
my $vepcache  = "/storage1/fs1/songcao/Active/Database/hg38_database/vep/v102";

# resources used by pipeline
my $GNOMAD_VCF          = "/storage1/fs1/songcao/Active/Database/tonlydb/af-only-gnomad.hg38.vcf.gz";
my $PANEL_OF_NORMALS_VCF= "/storage1/fs1/dinglab/Active/Projects/austins2/tools/somaticwrapper_tonly.v1.0/db/gatk/GDC-gatk4_panel_of_normals/6c4c4a48-3589-4fc0-b1fd-ce56e88c06e4/gatk4_mutect2_4136_pon.hg38.vcf.gz";
my $COMMON_BIALLELIC    = "/storage1/fs1/dinglab/Active/Projects/austins2/tools/somaticwrapper_tonly.v1.0/db/gatk/mutect2_gatk-best-practices.broadIns/af-only-gnomad.hg38.common_biallelic.chr1-22XY.vcf";
my $DB_SNP_NO_COSMIC    = "/storage1/fs1/songcao/Active/Database/hg38_database/cosmic/00-All.HG38.pass.cosmic.vcf";
my $f_ref_annot         = "/storage1/fs1/songcao/Active/Database/hg38_database/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa";

# infer chr prefix from ref (optional convenience)
my $first_line = `head -n 1 $h38_REF`;
if ($first_line =~ /^\>chr/) { $chr_status = 1; }

# input samples
opendir(my $DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir $DH;
closedir $DH;

# ---------------- SUBS ----------------

sub bsub_mutect2 {
  @J1_JIDS = (); $J2_JID = ''; $J3_JID = '';

  my $sample_full_path = shift;
  my $sample_name      = shift;

  my $out_mutect2="$sample_full_path/mutect2";	
  if(! -d $out_mutect2) { `mkdir -p $out_mutect2`; } 

  my @chrlist=(1..22,'X','Y');
  foreach my $chr (@chrlist) {
    my $chr1 = $chr_status ? "chr$chr" : $chr;
    my $current_job_file = "j1_mutect2_${sample_name}_${chr1}.sh";
    my $IN_bam_T = "$sample_full_path/$sample_name.remDup.bam";
    my $lsf_out  = "$lsf_file_dir/$current_job_file.out";
    my $lsf_err  = "$lsf_file_dir/$current_job_file.err";
    my $f_out_gz = "$out_mutect2/$chr1-f1r2.tar.gz"; 
    my $f_out_vcf= "$out_mutect2/$chr1-unfiltered.vcf";

    `rm -f $lsf_out $lsf_err`;

    open(my $FH, ">", "$job_files_dir/$current_job_file") or die $!;
    print $FH "#!/bin/bash\n";
    print $FH "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK Mutect2 -I $IN_bam_T -R $h38_REF -L $chr1 --germline-resource $GNOMAD_VCF -pon $PANEL_OF_NORMALS_VCF --f1r2-tar-gz $f_out_gz -O $f_out_vcf\n";
    close $FH;

    my $cmd = "bsub -q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a 'docker(scao/dailybox)' -o $lsf_out -e $lsf_err bash $job_files_dir/$current_job_file";
    my $out = `$cmd`;
    if ($out =~ /Job <(\d+)>/) { push @J1_JIDS, $1; }
  }
}

sub bsub_filter_mutect2 {
  my $sample_full_path = shift;
  my $sample_name      = shift;

  my $current_job_file = "j2_filter_mutect2_${sample_name}.sh";
  my $out_mutect2 = "$sample_full_path/mutect2";
  my $lsf_out = "$lsf_file_dir/$current_job_file.out";
  my $lsf_err = "$lsf_file_dir/$current_job_file.err";
  `rm -f $lsf_out $lsf_err`;

  my @chrlist=(1..22,'X','Y');
  my $all_unfiltered_input = join(" ", map {
    my $c = $chr_status ? "chr$_" : $_;
    "I=$out_mutect2/$c-unfiltered.vcf"
  } @chrlist);

  my $all_unfiltered_stats_input = join(" ", map {
    my $c = $chr_status ? "chr$_" : $_;
    "-stats $out_mutect2/$c-unfiltered.vcf.stats"
  } @chrlist);

  my $all_f1r2_input = join(" ", map {
    my $c = $chr_status ? "chr$_" : $_;
    "-I $out_mutect2/$c-f1r2.tar.gz"
  } @chrlist);

  my $IN_bam_T = "$sample_full_path/$sample_name.remDup.bam";
  my $f_merged_vcf   = "$out_mutect2/merged-unfiltered.vcf";
  my $f_merged_status= "$out_mutect2/merged-unfiltered.vcf.stats"; 
  my $f_ori          = "$out_mutect2/read-orientation-model.tar.gz";
  my $f_sum          = "$out_mutect2/getpileupsummaries.table"; 
  my $f_tab          = "$out_mutect2/contamination.table"; 
  my $f_tab2         = "$out_mutect2/contamination2.table";
  my $f_seg          = "$out_mutect2/segments.table";	
  my $f_filtered_vcf = "$out_mutect2/filtered.vcf"; 

  open(my $FH, ">", "$job_files_dir/$current_job_file") or die $!;
  print $FH "#!/bin/bash\n";
  print $FH "$java_bin -Xmx16g -jar $picardexe GatherVcfs $all_unfiltered_input O=$f_merged_vcf\n";
  print $FH "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK MergeMutectStats $all_unfiltered_stats_input -O $f_merged_status\n"; 
  print $FH "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK LearnReadOrientationModel $all_f1r2_input -O $f_ori\n";
  print $FH "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx8g -jar $GATK GetPileupSummaries -I $IN_bam_T -V $COMMON_BIALLELIC -L $COMMON_BIALLELIC -O $f_sum\n";
  print $FH "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK CalculateContamination -I $f_sum --tumor-segmentation $f_seg -O $f_tab\n";
  print $FH "    ".$run_script_path."check_contamination.pl $f_tab $f_tab2\n"; 
  print $FH "$java_bin -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -Xmx16g -jar $GATK FilterMutectCalls -V $f_merged_vcf -R $h38_REF --tumor-segmentation $f_seg --contamination-table $f_tab2 --ob-priors $f_ori -O $f_filtered_vcf\n";	
  close $FH;

  my $dep = @J1_JIDS ? "-w \"" . (join(" && ", map { "done($_)" } @J1_JIDS)) . "\" " : "";
  my $cmd = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub ${dep}-q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a 'docker(scao/dailybox)' -o $lsf_out -e $lsf_err bash $job_files_dir/$current_job_file";
  my $out = `$cmd`;
  if ($out =~ /Job <(\d+)>/) { $J2_JID = $1; }
}

sub bsub_parse_mutect2 {
  my $sample_full_path = shift;
  my $sample_name      = shift;

  my $current_job_file = "j3_parse_mutect2_${sample_name}.sh";
  my $out_mutect2 = "$sample_full_path/mutect2";
  my $f_filtered_vcf = "$out_mutect2/filtered.vcf";
  my $f_passed_vcf   = "$out_mutect2/filtered.pass.vcf";
  my $lsf_out = "$lsf_file_dir/$current_job_file.out";
  my $lsf_err = "$lsf_file_dir/$current_job_file.err";
  `rm -f $lsf_out $lsf_err`;
	
  open(my $FH, ">", "$job_files_dir/$current_job_file") or die $!;
  print $FH "#!/bin/bash\n";
  print $FH "export JAVA_HOME=$java_dir\n";
  print $FH "export JAVA_OPTS=\"-Xmx10g\"\n";
  print $FH "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
  print $FH "cat > $out_mutect2/mutect2_dbsnp_filter.input <<EOF\n";
  print $FH "mutect.dbsnp.annotator = $snpsift\n";
  print $FH "streka.dbsnp.db = $DB_SNP_NO_COSMIC\n";
  print $FH "streka.dbsnp.rawvcf = ./filtered.pass.vcf\n";
  print $FH "streka.dbsnp.mode = filter\n";
  print $FH "streka.dbsnp.passfile  = ./mutect2.somatic.dbsnp_pass.vcf\n";
  print $FH "streka.dbsnp.dbsnpfile = ./mutect2.somatic.dbsnp_present.vcf\n";
  print $FH "EOF\n";
  print $FH "cd $out_mutect2\n";
  print $FH "    ".$run_script_path."filter_mutect2.pl $f_filtered_vcf $f_passed_vcf $mincov_t $minvaf\n";
  print $FH "    ".$run_script_path."dbsnp_filter.pl ./mutect2_dbsnp_filter.input\n";
  close $FH;	

  my $dep2 = $J2_JID ? "-w \"done($J2_JID)\" " : "";
  my $cmd = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub ${dep2}-q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a 'docker(scao/dailybox)' -o $lsf_out -e $lsf_err bash $job_files_dir/$current_job_file";
  my $out = `$cmd`;
  if ($out =~ /Job <(\d+)>/) { $J3_JID = $1; }
}

sub bsub_vcf_2_maf {
  my $sample_full_path = shift;
  my $sample_name      = shift;

  my $current_job_file = "j4_vcf_2_maf.${sample_name}.sh";
  my $lsf_out = "$lsf_file_dir/$current_job_file.out";
  my $lsf_err = "$lsf_file_dir/$current_job_file.err";
  `rm -f $lsf_out $lsf_err`;

  my $f_m_vcf_in = "$sample_full_path/mutect2/mutect2.somatic.dbsnp_pass.vcf";

  open(my $FH, ">", "$job_files_dir/$current_job_file") or die $!;
  print $FH "#!/bin/bash\n";
  print $FH "F_VCF_1=$sample_full_path/merged.mutect2.vcf\n";
  print $FH "F_VCF_2=$sample_full_path/$sample_name.mutect2.vcf\n";
  print $FH "F_VEP_1=$sample_full_path/merged.VEP.mutect2.vcf\n";
  print $FH "F_VEP_2=$sample_full_path/$sample_name.mutect2.vep.vcf\n";
  print $FH "F_maf=$sample_full_path/$sample_name.mutect2.maf\n";
  print $FH "RUNDIR=$sample_full_path\n";
  print $FH "F_log=$sample_full_path/vep.merged.mutect2.log\n";
  print $FH "cat > \${RUNDIR}/vep.merged.mutect2.input <<EOF\n";
  print $FH "merged.vep.vcf = ./merged.mutect2.vcf\n";
  print $FH "merged.vep.output = ./merged.VEP.mutect2.vcf\n";
  print $FH "merged.vep.vep_cmd = $vepannot\n";
  print $FH "merged.vep.cachedir = $vepcache\n";
  print $FH "merged.vep.reffasta = $f_ref_annot\n";
  print $FH "merged.vep.assembly = GRCh38\n";
  print $FH "EOF\n";
  print $FH "rm -f \${F_log}\n";
  print $FH "cd $sample_full_path\n";
  print $FH "cp $f_m_vcf_in \${F_VCF_1}\n";  
  print $FH "    ".$run_script_path."vep_annotator.pl ./vep.merged.mutect2.input >&./vep.merged.mutect2.log\n";
  print $FH "rm -f \${F_VCF_2} \${F_VEP_2}\n";
  print $FH "ln -s \${F_VCF_1} \${F_VCF_2}\n";
  print $FH "ln -s \${F_VEP_1} \${F_VEP_2}\n";
  print $FH "    ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id ${sample_name}_T --normal-id ${sample_name}_N --ref-fasta $f_ref_annot --file-tsl /storage1/fs1/songcao/Active/Database/tsl/wgEncodeGencodeTranscriptionSupportLevelV23.txt\n";	
  close $FH;

  my $dep3 = $J3_JID ? "-w \"done($J3_JID)\" " : "";
  my $cmd = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub ${dep3}-g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a 'docker(ensemblorg/ensembl-vep:release_102.0)' -o $lsf_out -e $lsf_err bash $job_files_dir/$current_job_file";
  my $out = `$cmd`;
  if ($out =~ /Job <(\d+)>/) { push @ALL_J4_JIDS, $1; }
}

sub bsub_run_report {
  my $working_name = (split(/\//,$run_dir))[-1];
  my $current_job_file = "j5_Run_report_${working_name}.sh";
  my $lsf_out = "$lsf_file_dir/$current_job_file.out";
  my $lsf_err = "$lsf_file_dir/$current_job_file.err";
  `rm -f $lsf_out $lsf_err`;

  my $f_maf    = "$run_dir/$working_name.mutect2.maf";
  my $f_maf_rc = "$f_maf.rc";

  open(my $FH, ">", "$job_files_dir/$current_job_file") or die $!;
  print $FH "#!/bin/bash\n";
  print $FH "    ".$run_script_path."generate_final_report.pl $run_dir $status_exonic\n";
  print $FH "    ".$run_script_path."add_rc.pl $run_dir $f_maf $f_maf_rc\n";
  close $FH;

  my $dep_all = @ALL_J4_JIDS ? "-w \"" . (join(" && ", map { "done($_)" } @ALL_J4_JIDS)) . "\" " : "";
  my $cmd = "bsub ${dep_all}-q $q_name -g /$compute_username/$group_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a 'docker(scao/dailybox)' -o $lsf_out -e $lsf_err bash $job_files_dir/$current_job_file";
  my $out = `$cmd`;
  if ($out =~ /Job <(\d+)>/) { $J5_JID = $1; }
}

sub submit_all_for_sample {
  my ($sample_full_path, $sample_name) = @_;
  bsub_mutect2($sample_full_path, $sample_name);
  bsub_filter_mutect2($sample_full_path, $sample_name);
  bsub_parse_mutect2($sample_full_path, $sample_name);
  bsub_vcf_2_maf($sample_full_path, $sample_name);
}

# ---------------- MAIN ----------------

if ($step_number < 5) {
  for my $sample_name (@sample_dir_list) {
    next if ($sample_name =~ /\./ || $sample_name =~ /worklog/);
    my $sample_full_path = "$run_dir/$sample_name";
    next unless -d $sample_full_path;

    print $yellow, "\nSubmitting jobs for sample $sample_name ...",$normal,"\n";

    if    ($step_number == 0) { submit_all_for_sample($sample_full_path, $sample_name) }
    elsif ($step_number == 1) { bsub_mutect2($sample_full_path, $sample_name) }
    elsif ($step_number == 2) { bsub_filter_mutect2($sample_full_path, $sample_name) }
    elsif ($step_number == 3) { bsub_parse_mutect2($sample_full_path, $sample_name) }
    elsif ($step_number == 4) { bsub_vcf_2_maf($sample_full_path, $sample_name) }
  }

  if ($step_number == 0) {
    bsub_run_report();  # waits on ALL j4 jobs
    print "Queued run-level job: j5=$J5_JID\n";
  }
}

if ($step_number == 5) {
  print $yellow, "Submitting jobs for generating the run-level report ...",$normal,"\n";
  bsub_run_report();
}

exit 0;