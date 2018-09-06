# Perform vcf_filter.pl (VAF, read depth, and indel length) on strelka output
# Not performing dbSnP filtering
#
# Output is strelka/filter_out/strelka.somatic.snv.all.filtered.vcf

sub parse_strelka {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $filter_dir = shift;
    my $input_snv = shift;  
    my $strelka_vcf_filter_config = shift;

    die "Error: Input data file $input_snv does not exist\n" if (! -e $input_snv);

    my $filter_results = "$results_dir/strelka/filter_out";
    system("mkdir -p $filter_results");

# Create run script
    my $runfn = "$job_files_dir/j3_parse_strelka.sh";
    print STDERR "Writing to $runfn\n";
    open(OUT, ">$runfn") or die $!;

    my $vcf_filtered_fn = "$filter_results/strelka.somatic.snv.all.filtered.vcf";

# Run vcf_filter.py family of filters: VAF, read depth, and indel length
#    * Outputs strelka/filter_out/strelka.somatic.snv.all.filtered.vcf -> used for merge_vcf

    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx10g\"
export VARSCAN_DIR="/usr/local"

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
    my $cmd = "bash < $runfn\n";
    print STDERR "Executing:\n $cmd \n";

    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
