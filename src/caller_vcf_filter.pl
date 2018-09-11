# Perform post-call VCF filtering based on VAF, indel length, and read depth.
# Parses varscan, pindel, and strelka VCFs
# Filtering is performed by pyvcf-based scripts in ./vcf_filter
# Configuration file defines all adjustable filter parameters.
# Output filename is given by $output_vcf, written to $results_dir/caller_vcf_filter/$output_vcf

sub caller_vcf_filter {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $filter_dir = shift;
    my $input_vcf = shift;  
    my $output_vcf = shift;  #
    my $caller = shift;  # strelka, varscan, or pindel
    my $vcf_filter_config = shift;
    my $bypass = shift;  # boolean: will skip filtering if defined


    my $bypass_vcf = $bypass ? "--bypass" : "";
    die "Error: Input data file $input_vcf does not exist\n" if (! -e $input_vcf);
    die "Error: Caller not defined\n" if (! $caller);

    my $filter_results = "$results_dir/caller_vcf_filter";
    system("mkdir -p $filter_results");

# Create run script
    my $runfn = "$job_files_dir/j_caller_vcf_filter.sh";
    print STDERR "Writing to $runfn\n";
    open(OUT, ">$runfn") or die $!;

    my $vcf_filtered_fn = "$filter_results/$output_vcf";

# Run vcf_filter.py family of filters: VAF, read depth, and indel length

    print OUT <<"EOF";
#!/bin/bash

>&2 echo Running combined vcf_filter.py filters: VAF, read depth, and indel length
export PYTHONPATH="$filter_dir:\$PYTHONPATH"
bash $filter_dir/run_combined_vcf_filter.sh $dbsnp_filtered_fn $caller $vcf_filter_config $vcf_filtered_fn $bypass_vcf
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

EOF
    close OUT;
    my $cmd = "bash < $runfn\n";
    print STDERR "Executing:\n $cmd \n";

    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
