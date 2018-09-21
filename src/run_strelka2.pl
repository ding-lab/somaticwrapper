# Run Strelka

    # Strelka 2 results: $results_dir/strelka2/strelka_out/results/variants/somatic.snvs.vcf.gz

# Optional arg manta_vcf is explained in best practice here; https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#configuration
# implemented only for is_strelka2
sub run_strelka2 {
    my $IN_bam_T = shift;
    my $IN_bam_N = shift;
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $strelka2_dir = shift;  
    my $ref = shift;
    my $strelka_config = shift;
    my $manta_vcf = shift;    

    print STDERR "Running Strelka 2\n";
    $strelka_bin="$strelka2_dir/bin/configureStrelkaSomaticWorkflow.py";

    die "Error: Tumor BAM $IN_bam_T does not exist\n" if (! -e $IN_bam_T);
    die "Error: Tumor BAM $IN_bam_T is empty\n" if (! -s $IN_bam_T);
    die "Error: Normal BAM $IN_bam_N does not exist\n" if (! -e $IN_bam_N);
    die "Error: Normal BAM $IN_bam_N is empty\n" if (! -s $IN_bam_N);

    my $strelka_out=$results_dir."/strelka2/strelka_out";

    # Read configuration file into %params
    # Same format as used for varscan 
    my %params = get_config_params($strelka_config, 0);

    # currently strelka_flags used only for strelka2
    my $strelka_flags = "";
    if ($params{'is_exome'}) {
        $strelka_flags .= " --exome ";
    }
    if ($manta_vcf) {
        $strelka_flags .= " --indelCandidates $manta_vcf ";
    }

    my $expected_out;

    my $runfn = "$job_files_dir/j1_streka.sh"; 
    print STDERR "Writing to $runfn\n";
    open(OUT, ">$runfn") or die $!;

#
# strelka 2
#
    print STDERR "Executing Strelka 2\n";
    print OUT <<"EOF";
#!/bin/bash

if [ -d $strelka_out ] ; then
    rm -rf $strelka_out
fi

$strelka_bin $strelka_flags --normalBam $IN_bam_N --tumorBam $IN_bam_T --referenceFasta $ref --config $strelka_config --runDir $strelka_out
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


cd $strelka_out
ls
./runWorkflow.py -m local -j 8 -g 4
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

EOF
    close OUT;

    $expected_out="$strelka_out/results/variants/somatic.snvs.vcf.gz";

    my $cmd = "bash < $runfn\n";

    print STDERR $cmd."\n";
    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;

    printf STDERR "Testing output $expected_out\n";
    die "Error: Did not find expected output file $expected_out\n" if (! -e $expected_out);
    printf STDERR "OK\n";
}

1;
