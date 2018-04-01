# bam readcount not used because not using FP Filter
#my $bam_readcount = "/usr/local/bin/bam-readcount";

# Pre-CWL arrangement:
    # The following files are created by prior steps, $sample_full_path/strelka/strelka_out/results and accessed here
        # passed.somatic.indels.vcf -> renamed to strelka.somatic.indel.strlk_pass.gvip.vcf with genomevip_label
        # passed.somatic.snvs.vcf   -> renamed to strelka.somatic.snv.strlk_pass.gvip.vcf with genomevip_label

    # Output of this script will be in strelka/filter_out
    # * strelka.somatic.indel.strlk_pass.vcf  (vep annotation)
    # * strelka.somatic.snv.strlk_pass.vcf    (vep annotation)
    # * strelka.somatic.snv.all.dbsnp_pass.vcf (vep annotation and merge_vcf)

# With CWL, we need to pass files we are operating on explicitly as input argument
#  * Only passed.somatic.snvs.vcf is actually used.  We will keep this as $input_snv
#  * Skipping genomevip_label annotation.  Not dealing with passed.somatic.indels.vcf here at all
#  * only output is strelka.somatic.snv.all.dbsnp_pass.vcf 

sub parse_strelka {
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $input_snv = shift;  # New to CWL: pass this filename explicitly (passed.somatic.snvs.vcf)

    die "Error: dbSnP database file $dbsnp_db does not exist\n" if (! -e $dbsnp_db);

    $current_job_file = "j3_parse_strelka.sh";

    my $bsub = "bash";
    my $strelka_results = "$sample_full_path/strelka/strelka_out/results";
    my $filter_results = "$sample_full_path/strelka/filter_out";
    system("mkdir -p $filter_results");

# create strelka_dbsnp_filter.snv.input
    my $dbsnp_config = "$filter_results/strelka_dbsnp_filter.snv.input";
    print("Writing to $dbsnp_config\n");
    open(OUT, ">$dbsnp_config") or die $!;
    print OUT <<"EOF";
streka.dbsnp.snv.annotator = $snpsift_jar
streka.dbsnp.snv.db = $dbsnp_db
streka.dbsnp.snv.rawvcf = $input_snv
streka.dbsnp.snv.mode = filter
streka.dbsnp.snv.passfile  = $filter_results/strelka.somatic.snv.all.dbsnp_pass.vcf
streka.dbsnp.snv.dbsnpfile = $filter_results/strelka.somatic.snv.all.dbsnp_present.vcf
EOF

# Create run script
    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;


# Note that dbsnp_filter.pl automatically adds dbsnp_anno.vcf suffix to rawvcf when creating output
# Step 5 (dbSnP Filter on SNV) creates these two files:
#   strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf  -> used for merge_vcf
#   strelka/filter_out/strelka.somatic.snv.all.dbsnp_present.vcf 

    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx10g\"
export VARSCAN_DIR="/usr/local"

# step 5
$perl $gvip_dir/dbsnp_filter.pl $dbsnp_config

EOF
    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
