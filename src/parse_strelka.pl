my $datd="/data/B_Filter";
my $jar="/usr/local/snpEff/SnpSift.jar";

# filtered database created in B_Filter
my $db="$datd/dbsnp.noCOSMIC.vcf.gz";
#my $db="$datd/short.dbsnp.noCOSMIC.vcf.gz";

# The following files are created by prior steps, $sample_full_path/strelka/strelka_out/results
    # all.somatic.indels.vcf
    # all.somatic.snvs.vcf
    # passed.somatic.indels.vcf
    # passed.somatic.snvs.vcf
    # strelka_dbsnp_filter.indel.input
    # strelka_dbsnp_filter.snv.input
    # strelka_fpfilter.snv.input

# Output of this script will be in $sample_full_path/filter_out


sub parse_strelka{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;

    $current_job_file = "j3_parse_strelka".$sample_name.".sh";

    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
#    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $strelka_results = "$sample_full_path/strelka/strelka_out/results";
    my $filter_results = "$sample_full_path/strelka/filter_out";
    system("mkdir -p $filter_results");

# create strelka_dbsnp_filter.snv.input
    my $dbsnp_snv = "$filter_results/strelka_dbsnp_filter.snv.input";
    print("Writing to $dbsnp_snv\n");
    open(OUT, ">$dbsnp_snv") or die $!;
    print OUT <<"EOF";
streka.dbsnp.snv.annotator = $jar
streka.dbsnp.snv.db = $db
streka.dbsnp.snv.rawvcf = $filter_results/strelka.somatic.snv.strlk_pass.gvip.vcf
streka.dbsnp.snv.mode = filter
streka.dbsnp.snv.passfile  = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
streka.dbsnp.snv.dbsnpfile = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_present.vcf
EOF


# create strelka_dbsnp_filter.indel.input
    my $dbsnp_indel = "$filter_results/strelka_dbsnp_filter.indel.input";
    print("Writing to $dbsnp_indel\n");
    open(OUT, ">$dbsnp_indel") or die $!;
    print OUT <<"EOF";
streka.dbsnp.indel.annotator = $jar
streka.dbsnp.indel.db = $db
streka.dbsnp.indel.rawvcf = $filter_results/strelka.somatic.indel.strlk_pass.gvip.vcf
streka.dbsnp.indel.mode = filter
streka.dbsnp.indel.passfile  = $filter_results/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf
streka.dbsnp.indel.dbsnpfile = $filter_results/strelka.somatic.indel.all.gvip.dbsnp_present.vcf
EOF


# create strelka_fpfilter.snv.input
    my $fp_snv = "$filter_results/strelka_fpfilter.snv.input";
    print("Writing to $fp_snv\n");
    open(OUT, ">$fp_snv") or die $!;
    print OUT <<"EOF";
strelka.fpfilter.snv.bam_readcount = /usr/local/bin/bam-readcount
strelka.fpfilter.snv.bam_file = $IN_bam_T
strelka.fpfilter.snv.REF = $REF
strelka.fpfilter.snv.variants_file = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
strelka.fpfilter.snv.passfile = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_pass.vcf
strelka.fpfilter.snv.failfile = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_fail.vcf
strelka.fpfilter.snv.rc_in = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.in.vcf
strelka.fpfilter.snv.rc_out = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.out.vcf
strelka.fpfilter.snv.fp_out = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp.out.vcf
strelka.fpfilter.snv.min_mapping_qual = 0
strelka.fpfilter.snv.min_base_qual = 15
strelka.fpfilter.snv.min_num_var_supporting_reads = 4
strelka.fpfilter.snv.min_var_allele_freq = 0.05
strelka.fpfilter.snv.min_avg_rel_read_position = 0.10
strelka.fpfilter.snv.min_avg_rel_dist_to_3prime_end = 0.10
strelka.fpfilter.snv.min_var_strandedness = 0.01
strelka.fpfilter.snv.min_allele_depth_for_testing_strandedness = 5
strelka.fpfilter.snv.min_ref_allele_avg_base_qual = 30
strelka.fpfilter.snv.min_var_allele_avg_base_qual = 30
strelka.fpfilter.snv.max_rel_read_length_difference = 0.25
strelka.fpfilter.snv.max_mismatch_qual_sum_for_var_reads = 150
strelka.fpfilter.snv.max_avg_mismatch_qual_sum_difference = 150
strelka.fpfilter.snv.min_ref_allele_avg_mapping_qual = 30
strelka.fpfilter.snv.min_var_allele_avg_mapping_qual = 30
strelka.fpfilter.snv.max_avg_mapping_qual_difference = 50
EOF

# Create run script
    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;


# Note that dbsnp_filter.pl automatically adds dbsnp_anno.vcf suffix to rawvcf when creating output
# Step 5 creates these two files:
#   strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf  - not empty
#   strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_present.vcf - empty (header only)
# Step 6 creates these two files:
#   strelka/filter_out/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf  - not empty
#   strelka/filter_out/strelka.somatic.indel.all.gvip.dbsnp_present.vcf - empty (header only)
# Step 7 creates this file:
#   strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.fp.out.vcf - not empty

    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx512m\"
export VARSCAN_DIR="/usr/local"

$perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/all.somatic.snvs.vcf $filter_results/strelka.somatic.snv.all.gvip.vcf
$perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/all.somatic.indels.vcf $filter_results/strelka.somatic.indel.all.gvip.vcf
$perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/passed.somatic.snvs.vcf $filter_results/strelka.somatic.snv.strlk_pass.gvip.vcf
$perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/passed.somatic.indels.vcf $filter_results/strelka.somatic.indel.strlk_pass.gvip.vcf
$perl $gvip_dir/dbsnp_filter.pl $filter_results/strelka_dbsnp_filter.snv.input
$perl $gvip_dir/dbsnp_filter.pl $filter_results/strelka_dbsnp_filter.indel.input
$perl $gvip_dir/snv_filter.pl $filter_results/strelka_fpfilter.snv.input

EOF
    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");
    
    system ( $bsub_com ); 

}

1;
