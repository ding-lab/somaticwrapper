# Apply dbSnP filter to a VCF
#
# dbSnP filter processes a VCF and removes/retains variants based on whether they are present
# in a variant database (e.g., dbSnP)
# Given input data "/path/to/data/INPUT.vcf",
#   * Create '$results_dir/dbsnp_filter/INPUT.dbsnp_pass.vcf' with all entries which are NOT present in variant db
#   * Create '$results_dir/dbsnp_filter/INPUT.dbsnp_present.vcf' with all entries which ARE present in variant db
# If $bypass is true, then dbsnp_pass.vcf will have ALL variants.  Note that VCF will be annotated with variant names
#   (ID column) obtained from dbsnp_db
# If dbsnp_db is not defined, skip processing altogether and the output passfile is simply a link to the input

use File::Basename;

sub dbsnp_filter {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $input_vcf = shift;
    my $bypass = shift;

    
    my $filter_results = "$results_dir/dbsnp_filter";
    print STDERR "Filter results: $filter_results\n";
    system("mkdir -p $filter_results");

    # Strip the path and extension from $input_vcf, and create the dbsnp_pass and dbsnp_present files from it
    my($filename, $dirs, $suffix) = fileparse($input_vcf);
    my $pass_fn = $filter_results . "/" . $filename . ".dbsnp_pass.vcf";
    my $present_fn = $filter_results . "/" . $filename . ".dbsnp_present.vcf";

    # Test if dbsnp_db is defined and whether it is a file
    if ($dbsnp_db eq "") {
        # If dbsnp_db not defined, simply link output to input and return
        my $infn = `readlink -f $input_vcf`;
        chomp $infn;
        print STDERR "Skipping dbSnP filter because dbsnp_db not defined.  Creating $pass_fn as link to $infn\n";
        my $return_code = system ( "ln -fs $infn $pass_fn" );
        die("Exiting ($return_code).\n") if $return_code != 0;
        return;
    } else {
        die "Error: dbSnP database file $dbsnp_db does not exist\n" if (! -e $dbsnp_db);
    }

    my $mode_str = "filter";
    if ($bypass) {  # this is similar to not defining dbsnp_db in that there is no filtering, but dbsnp_filter.pl is called, and output VCF is annotated with variant names
        $mode_str = "annotate";
    }
    # $dbsnp_filtered_snv_fn is the output of dbsnp filter of SNV calls
    my $config_fn = "$filter_results/dbsnp_filter.config";
    print STDERR "Writing to $config_fn\n";
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.snv.annotator = $snpsift_jar
varscan.dbsnp.snv.db = $dbsnp_db
varscan.dbsnp.snv.rawvcf = $input_vcf
varscan.dbsnp.snv.mode = $mode_str
varscan.dbsnp.snv.passfile  = $pass_fn
varscan.dbsnp.snv.presentfile = $present_fn
EOF

    my $runfn = "$job_files_dir/j_dbsnp_filter.sh";
    print STDERR "Writing to $runfn\n";
    open(OUT, ">$runfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xms256m -Xmx10g\"

$perl $gvip_dir/dbsnp_filter.pl  $config_fn
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
