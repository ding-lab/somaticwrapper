# Process pindel run output and generate VCF.

# Calls GenomeVIP/pindel_filter.pl, which performs the following:
# 1. apply CvgVafStrand Filter (coverage) to pindel output
# 2. Convert reads to VCF
# 3. apply homopolymer filter
# Note that caller vcf filtering (length, depth, VAF) is not performed by this step
# 
# Output:
# $results_dir/pindel/filter_out/$pindel_raw.CvgVafStrand_pass.Homopolymer_pass.vcf
#
# if bypass is true, skip filtering for CvgVafStrand and Homopolymer in pindel_filter, and just output VCF file

sub parse_pindel {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $reference = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $pindel_dir = shift;
    my $pindel_config = shift;
    my $pindel_raw_in = shift; 
    my $no_delete_temp = shift;
    my $bypass = shift;

    if (! $no_delete_temp) {
        $no_delete_temp = 0; 
    }

    my $filter_results = "$results_dir/pindel/filter_out";
    system("mkdir -p $filter_results");

    # pindel_filter writes all output data to the same directory as input data.
    # To get around this, make link to input data in output directory.  make_data_link()
    # defined in parse_varscan_snv.pl
    my $pindel_raw = make_data_link($pindel_raw_in, $filter_results);

    # bypass flag will skip two filters in parse_pindel, and vcf_filter.  
    my $bypass_str = $bypass ? "pindel.filter.skip_filter1 = true\npindel.filter.skip_filter2 = true" : "";

    # Set up parse_pindel configuration file.  It takes $pindel_config file and adds several lines to it
    die "$pindel_config does not exist\n" unless (-f $pindel_config);
    my $config_fn = "$filter_results/pindel_filter.input";
    print STDERR "Copying $pindel_config to $config_fn and appending\n";
    $errcode = system("cp $pindel_config $config_fn");
    die ("Error executing: $cmd \n $! \n") if ($errcode);
    open(OUT, ">>$config_fn") or die $!;
    print OUT <<"EOF";
pindel.filter.pindel2vcf = $pindel_dir/pindel2vcf
pindel.filter.variants_file = $pindel_raw
pindel.filter.REF = $reference
pindel.filter.date = 000000
$bypass_str
EOF

    # Run pindel_filter produces:
    #    pindel.out.raw.CvgVafStrand_pass 
    #    pindel.out.raw.CvgVafStrand_fail
    #    pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf  -> this is input into dbSnP filter
    #    pindel.out.raw.CvgVafStrand_pass.Homopolymer_fail.vcf  
    my $pindel_filter_out="$pindel_raw.CvgVafStrand_pass.Homopolymer_pass.vcf";  # This is the output of the script
    my $pindel_filter_cmd = "$perl $gvip_dir/pindel_filter.pl $config_fn";

    my $outfn = "$job_files_dir/j7_parse_pindel.sh";
    print STDERR "Writing to $outfn\n";
    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash

>&2 echo Running pindel_filter.pl
$pindel_filter_cmd
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

# Optionally delete intermediate files
#    - specifically, files with "_fail" in the filename
if [[ $no_delete_temp == 1 ]]; then
    >&2 echo Not deleting intermediate files
else
    >&2 echo Deleting intermediate \\"filter fail\\" files
    cd $filter_results
    rm -f \*_fail\* 
    tmp_base=\$(basename \$TMP)
    rm -f \$tmp_base
fi

EOF
    close OUT;
    my $cmd = "bash < $outfn\n";
    print STDERR "Executing:\n $cmd \n";

    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
