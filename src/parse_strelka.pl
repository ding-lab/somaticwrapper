
# Pass files we are operating on explicitly as input argument
#  * Only passed.somatic.snvs.vcf is actually used.  We will keep this as $input_snv
#  * Skipping genomevip_label annotation.  Not dealing with passed.somatic.indels.vcf here at all
#  * only output is strelka.somatic.snv.all.dbsnp_pass.filtered.vcf 

sub parse_strelka {
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $filter_dir = shift;
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $input_snv = shift;  # New to CWL: pass this filename explicitly (passed.somatic.snvs.vcf)
    my $strelka_vcf_filter_config = shift;

    # It would be helpul to allow dbsnp_db to be not set, which would imply skipping the filtering step
    # This is currently not supported: require dbsnp_db to be defined and a file
    if ($dbsnp_db eq "") {
        die("Error: dbsnp_db not defined\n");
    } else {
        die "Error: dbSnP database file $dbsnp_db does not exist\n" if (! -e $dbsnp_db);
    }
    die "Error: Input data file $input_snv does not exist\n" if (! -e $input_snv);

    $current_job_file = "j3_parse_strelka.sh";

    my $bsub = "bash";
    my $strelka_results = "$sample_full_path/strelka/strelka_out/results";
    my $filter_results = "$sample_full_path/strelka/filter_out";
    system("mkdir -p $filter_results");

# create strelka_dbsnp_filter.snv.input
    my $dbsnp_filtered_fn = "$filter_results/strelka.somatic.snv.all.dbsnp_pass.vcf";
    my $dbsnp_config = "$filter_results/strelka_dbsnp_filter.snv.input";
    print STDERR "Writing to $dbsnp_config\n";
    open(OUT, ">$dbsnp_config") or die $!;
    print OUT <<"EOF";
streka.dbsnp.snv.annotator = $snpsift_jar
streka.dbsnp.snv.db = $dbsnp_db
streka.dbsnp.snv.rawvcf = $input_snv
streka.dbsnp.snv.mode = filter
streka.dbsnp.snv.passfile  = $dbsnp_filtered_fn
streka.dbsnp.snv.dbsnpfile = $filter_results/strelka.somatic.snv.all.dbsnp_present.vcf
EOF

# Create run script
    my $outfn = "$job_files_dir/$current_job_file";
    print STDERR "Writing to $outfn\n";
    open(OUT, ">$outfn") or die $!;

    my $vcf_filtered_fn = "$filter_results/strelka.somatic.snv.all.dbsnp_pass.filtered.vcf";

# Note that dbsnp_filter.pl automatically adds dbsnp_anno.vcf suffix to rawvcf when creating output
# Step 5 (dbSnP Filter on SNV) creates these two files:
#   strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf  
#   strelka/filter_out/strelka.somatic.snv.all.dbsnp_present.vcf 
# Run vcf_filter.py family of filters: VAF, read depth, and indel length
#    * Reads strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf
#    * Outputs strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.filtered.vcf -> used for merge_vcf

    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx10g\"
export VARSCAN_DIR="/usr/local"

# dnsnp filtering
$perl $gvip_dir/dbsnp_filter.pl $dbsnp_config
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


>&2 echo Running combined vcf_filter.py filters: VAF, read depth, and indel length
export PYTHONPATH="$filter_dir:\$PYTHONPATH"
bash $filter_dir/run_combined_vcf_filter.sh $dbsnp_filtered_fn strelka $strelka_vcf_filter_config $vcf_filtered_fn
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


EOF
    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print STDERR "Executing:\n $bsub_com \n";

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
