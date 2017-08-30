
sub bsub_parse_strelka{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $script_dir = shift;

    $current_job_file = "j3_parse_strelka".$sample_name.".sh";

    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $strelka_results = "$sample_full_path/strelka/strelka_out/results";
    my $strelka_out = $sample_full_path."/strelka";

# create strelka_dbsnp_filter.snv.input
    my $dbsnp_snv = "$strelka_out/strelka_out/results/strelka_dbsnp_filter.snv.input";
    print("Writing to $dbsnp_snv\n");
    open(OUT, ">$dbsnp_snv") or die $!;
    print OUT <<"EOF";
streka.dbsnp.snv.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar
streka.dbsnp.snv.db = /gscmnt/gc3027/dinglab/medseq/cosmic/00-All.brief.pass.cosmic.vcf
streka.dbsnp.snv.rawvcf = ./strelka.somatic.snv.strlk_pass.gvip.vcf
streka.dbsnp.snv.mode = filter
streka.dbsnp.snv.passfile  = ./strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
streka.dbsnp.snv.dbsnpfile = ./strelka.somatic.snv.all.gvip.dbsnp_present.vcf
EOF


# create strelka_dbsnp_filter.indel.input
    my $dbsnp_indel = "$strelka_out/strelka_out/results/strelka_dbsnp_filter.indel.input";
    print("Writing to $dbsnp_indel\n");
    open(OUT, ">$dbsnp_indel") or die $!;
    print OUT <<"EOF";
streka.dbsnp.indel.annotator = /gscmnt/gc2525/dinglab/rmashl/Software/bin/snpEff/20150522/SnpSift.jar
streka.dbsnp.indel.db = /gscmnt/gc3027/dinglab/medseq/cosmic/00-All.brief.pass.cosmic.vcf
streka.dbsnp.indel.rawvcf = ./strelka.somatic.indel.strlk_pass.gvip.vcf
streka.dbsnp.indel.mode = filter
streka.dbsnp.indel.passfile  = ./strelka.somatic.indel.all.gvip.dbsnp_pass.vcf
streka.dbsnp.indel.dbsnpfile = ./strelka.somatic.indel.all.gvip.dbsnp_present.vcf
EOF


# create strelka_fpfilter.snv.input
    my $fp_snv = "$strelka_out/strelka_out/results/strelka_fpfilter.snv.input";
    print("Writing to $fp_snv\n");
    open(OUT, ">$fp_snv") or die $!;
    print OUT <<"EOF";
strelka.fpfilter.snv.bam_readcount = /gscmnt/gc2525/dinglab/rmashl/Software/bin/bam-readcount/0.7.4/bam-readcount
strelka.fpfilter.snv.bam_file = $IN_bam_T
strelka.fpfilter.snv.REF = $REF
strelka.fpfilter.snv.variants_file = $sample_full_path/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
strelka.fpfilter.snv.passfile = $sample_full_path/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_pass.vcf
strelka.fpfilter.snv.failfile = $sample_full_path/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp_fail.vcf
strelka.fpfilter.snv.rc_in = $sample_full_path/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.in.vcf
strelka.fpfilter.snv.rc_out = $sample_full_path/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.rc.out.vcf
strelka.fpfilter.snv.fp_out = $sample_full_path/strelka/strelka_out/results/strelka.somatic.snv.all.gvip.dbsnp_pass.fp.out.vcf
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


    print OUT <<"EOF";
#!/bin/bash

cd $strelka_results

$perl $script_dir/genomevip_label.pl Strelka ./all.somatic.snvs.vcf ./strelka.somatic.snv.all.gvip.vcf
$perl $script_dir/genomevip_label.pl Strelka ./all.somatic.indels.vcf ./strelka.somatic.indel.all.gvip.vcf
$perl $script_dir/genomevip_label.pl Strelka ./passed.somatic.snvs.vcf ./strelka.somatic.snv.strlk_pass.gvip.vcf
$perl $script_dir/genomevip_label.pl Strelka ./passed.somatic.indels.vcf ./strelka.somatic.indel.strlk_pass.gvip.vcf
$perl $script_dir/dbsnp_filter.pl ./strelka_dbsnp_filter.snv.input
$perl $script_dir/dbsnp_filter.pl ./strelka_dbsnp_filter.indel.input
$perl $script_dir/snv_filter.pl ./strelka_fpfilter.snv.input

EOF
    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");
    
    print("Aborting\n");
die();

    system ( $bsub_com ); 

}

1;
