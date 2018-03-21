# bam readcount not used because not using FP Filter
#my $bam_readcount = "/usr/local/bin/bam-readcount";

# The following files are created by prior steps, $sample_full_path/strelka/strelka_out/results
    # all.somatic.indels.vcf
    # all.somatic.snvs.vcf
    # passed.somatic.indels.vcf
    # passed.somatic.snvs.vcf

# Output of this script will be in strelka/filter_out


sub parse_strelka {
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $ref_dict = shift;  # this ends up being not used
    my $perl = shift;
    my $gvip_dir = shift;
    my $db = shift;
    my $snpsift_jar = shift;

    $current_job_file = "j3_parse_strelka.sh";

    my $bsub = "bash";
    my $strelka_results = "$sample_full_path/strelka/strelka_out/results";
    my $filter_results = "$sample_full_path/strelka/filter_out";
    system("mkdir -p $filter_results");

# create strelka_dbsnp_filter.snv.input
    my $dbsnp_snv = "$filter_results/strelka_dbsnp_filter.snv.input";
    print("Writing to $dbsnp_snv\n");
    open(OUT, ">$dbsnp_snv") or die $!;
    print OUT <<"EOF";
streka.dbsnp.snv.annotator = $snpsift_jar
streka.dbsnp.snv.db = $db
streka.dbsnp.snv.rawvcf = $filter_results/strelka.somatic.snv.strlk_pass.gvip.vcf
streka.dbsnp.snv.mode = filter
streka.dbsnp.snv.passfile  = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
streka.dbsnp.snv.dbsnpfile = $filter_results/strelka.somatic.snv.all.gvip.dbsnp_present.vcf
EOF

# Create run script
    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;


# Note that dbsnp_filter.pl automatically adds dbsnp_anno.vcf suffix to rawvcf when creating output
# Step 5 (dbSnP Filter on SNV) creates these two files:
#   strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf  -> used for merge_vcf
#   strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_present.vcf 
# Step 6 (dbSnP Filter on Indel) creates these two files - not used:
#   strelka/filter_out/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf  
#   strelka/filter_out/strelka.somatic.indel.all.gvip.dbsnp_present.vcf 
# Step 7 (FP Filter on SNV) creates this file - not used:
#   strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.fp.out.vcf 

# Note that in subsequent steps (merge_vcf) only the file strelka.somatic.snv.all.gvip.dbsnp_pass.vcf is used.
# output of steps 6 and 7 is discarded.  We are removing these steps at Song's suggestion
# We subsequently sort the VCF to make into a standard format for GATK merging

# TODO:  run genomevip_label at end of run steps, rather than beginning of parse steps.  This will make it easier
# to understand the relationship between the steps.

    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx10g\"
export VARSCAN_DIR="/usr/local"

# steps 1-4.  Only step 3 used for merging, although run_vep also uses 4.
# $perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/all.somatic.snvs.vcf $filter_results/strelka.somatic.snv.all.gvip.vcf
# $perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/all.somatic.indels.vcf $filter_results/strelka.somatic.indel.all.gvip.vcf
$perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/passed.somatic.snvs.vcf $filter_results/strelka.somatic.snv.strlk_pass.gvip.vcf
$perl $gvip_dir/genomevip_label.pl Strelka $strelka_results/passed.somatic.indels.vcf $filter_results/strelka.somatic.indel.strlk_pass.gvip.vcf

# step 5
$perl $gvip_dir/dbsnp_filter.pl $filter_results/strelka_dbsnp_filter.snv.input

# Following steps 6-7 are not used
#$perl $gvip_dir/dbsnp_filter.pl $filter_results/strelka_dbsnp_filter.indel.input
#$perl $gvip_dir/snv_filter.pl $filter_results/strelka_fpfilter.snv.input

EOF
    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 

}

1;
