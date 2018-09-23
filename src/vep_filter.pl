# Apply AF (allele frequency) and classification (e.g., exon-only) filter.
#
#   --bypass_af will skip AF filter by retaining all reads
#   --bypass_classification will skip classification filter by retaining all reads
#   --bypass will skip all reads

# Output is $results_dir/vep/vep_filtered.vcf

sub vep_filter {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $input_vcf = shift;  # Name of input VCF to process
    my $af_filter_config = shift;
    my $classification_filter_config = shift;
    my $bypass_af = shift;
    my $bypass_classification = shift;
    my $bypass = shift;
    my $debug = shift;

    my $filter_results = "$results_dir/vep_filter";
    system("mkdir -p $filter_results");

    my $bypass_str = $bypass ? "--bypass" : "";
    $bypass_str = $bypass_af ? "--bypass_af $bypass_str" : "$bypass_str";
    $bypass_str = $bypass_classification ? "--bypass_classification $bypass_str" : "$bypass_str";
    my $debug_str = $debug ? "--debug" : "";

    # Check to make sure filter config files exist
    die "AF filter config $af_filter_config does not exist\n" unless (-e $af_filter_config);
    die "Classification filter config $classification_filter_config does not exist\n" unless (-e $classification_filter_config);

    my $filtered_output_fn = "$filter_results/vep_filtered.vcf";

    my $runfn = "$job_files_dir/j_vep_filter.sh";
    print STDERR "Writing to $runfn\n";
    open(OUT, ">$runfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx512m\"

>&2 echo Filtering by AF and classification
export PYTHONPATH="$SWpaths::filter_dir:\$PYTHONPATH"

bash $SWpaths::filter_dir/run_vep_filters.sh $input_vcf $af_filter_config $classification_filter_config $filtered_output_fn $bypass_str $debug_str

EOF

    close OUT;
    my $cmd = "bash < $runfn";
    print STDERR "Executing:\n $cmd \n";

    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;

    print STDERR "Final results written to $filtered_output_fn\n";
}

1;
