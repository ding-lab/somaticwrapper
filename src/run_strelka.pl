# Run Strelka 1

    # Strelka 1 results: $results_dir/strelka/strelka_out/results/passed.somatic.snvs.vcf

sub run_strelka {
    my $IN_bam_T = shift;
    my $IN_bam_N = shift;
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $ref = shift;
    my $strelka_config = shift;

    print STDERR "Running Strelka 1\n";
    $strelka_bin="$SWpaths::strelka_dir/bin/configureStrelkaWorkflow.pl";

    die "Error: Tumor BAM $IN_bam_T does not exist\n" if (! -e $IN_bam_T);
    die "Error: Tumor BAM $IN_bam_T is empty\n" if (! -s $IN_bam_T);
    die "Error: Normal BAM $IN_bam_N does not exist\n" if (! -e $IN_bam_N);
    die "Error: Normal BAM $IN_bam_N is empty\n" if (! -s $IN_bam_N);

    my $strelka_out=$results_dir."/strelka/strelka_out";

    # Read configuration file into %params
    # Same format as used for varscan 
    my %params = get_config_params($strelka_config, 0);

    # currently strelka_flags used only for strelka2
    my $strelka_flags = "";
    if ($params{'is_exome'}) {
        $strelka_flags .= " --exome ";
    }

    my $expected_out;

    my $runfn = "$job_files_dir/j1_streka.sh"; 
    print STDERR "Writing to $runfn\n";
    open(OUT, ">$runfn") or die $!;

#
# Strelka 1
#
    print STDERR "Executing Strelka 1\n";
    print OUT <<"EOF";
#!/bin/bash

if [ -d $strelka_out ] ; then
    rm -rf $strelka_out
fi

$strelka_bin --normal $IN_bam_N --tumor $IN_bam_T --ref $ref --config $strelka_config --output-dir $strelka_out
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


cd $strelka_out
make -j 16
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

EOF
    close OUT;
    $expected_out="$strelka_out/results/passed.somatic.snvs.vcf";

    my $cmd = "bash < $runfn\n";

    print STDERR $cmd."\n";
    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;

    printf STDERR "Testing output $expected_out\n";
    die "Error: Did not find expected output file $expected_out\n" if (! -e $expected_out);
    printf STDERR "OK\n";
}

1;
