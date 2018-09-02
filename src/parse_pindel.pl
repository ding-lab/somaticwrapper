
# principal output used in merging: pindel/filter_out/pindel.out.current_final.dbsnp_pass.filtered.vcf

# if apply_filter is 0, skip filtering for CvgVafStrand and Homopolymer in pindel_filter, and just output VCF file

# TODO: add pindel.filter.apply_filter = true to pindel_filter input
# TODO: make dbsnp_db optional; if does not exist, skip this filtering

# CWL changes:
# * genomevip_labeling removed
# * Unnecessary copy operation removed (pindel.out.current_final.Somatic.vcf)
# * All input filenames explicitly passed
# * The `grep ChrID` command is moved to the run_pindel step.  This changes the input into this script,
#   which is `pindel_raw`
# * Optionally delete files pindel-raw.CvgVafStrand_fail and pindel-raw.CvgVafStrand_pass.Homopolymer_fail.vcf

sub parse_pindel {
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $filter_dir = shift;
    my $pindel_dir = shift;
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $pindel_config = shift;
    my $pindel_raw_in = shift; # NEW
    my $no_delete_temp = shift;
    my $pindel_vcf_filter_config = shift;
    my $bypass = shift;

    if (! $no_delete_temp) {
        $no_delete_temp = 0; # avoid empty variables
    }

    $current_job_file = "j7_parse_pindel.sh";

    my $bsub = "bash";
    my $filter_results = "$sample_full_path/pindel/filter_out";
    system("mkdir -p $filter_results");

    die "Error: dbSnP database file $dbsnp_db does not exist\n" if (! -e $dbsnp_db);

    # pindel_filter is unfortunately designed such that all output data is written to the same directory as input data, and
    # no way is provided to change that.  Since input data is passed, and we need
    # to be able to control where data is written to, we must create a soft-link to input data in output 
    # directory.  Link has to be created with an absolute filename
    die "Error: Pindel raw input file $pindel_raw_in does not exist\n" if (! -e $pindel_raw_in);
    $pindel_raw_in = `readlink -f $pindel_raw_in`;
    chomp $pindel_raw_in;

    my $errcode = system ("ln -fs $pindel_raw_in $filter_results "); 
    die ("Error executing: $cmd \n $! \n") if ($errcode);

    my $pindel_raw=$filter_results . "/" . basename($pindel_raw_in) ;
    my $filter_out="$pindel_raw.CvgVafStrand_pass.Homopolymer_pass.vcf";

    # This is the principal result of pindel_filter
    my $bypass_str;
    if ($bypass) {
        $bypass_str = "pindel.filter.skip_filter1 = true\npindel.filter.skip_filter2 = true";
    } else {
        $bypass_str = "";
    }

## Pindel Filter - below is input into pindel_filter.pl
# lines below are added to data from $pindel_config
    die "$pindel_config does not exist\n" unless (-f $pindel_config);

    my $out = "$filter_results/pindel_filter.input";
    print STDERR "Copying $pindel_config to $out and appending\n";
    $errcode = system("cp $pindel_config $out");
    die ("Error executing: $cmd \n $! \n") if ($errcode);

    open(OUT, ">>$out") or die $!;
    print OUT <<"EOF";
pindel.filter.pindel2vcf = $pindel_dir/pindel2vcf
pindel.filter.variants_file = $pindel_raw
pindel.filter.REF = $REF
pindel.filter.date = 000000
$bypass_str
EOF

## dbSnP Filter
# TODO: skip this step if $dbsnp_db not defined.  Not yet implemented
    my $dbsnp_filtered_fn = "$filter_results/pindel.out.current_final.dbsnp_pass.vcf";
    my $out = "$filter_results/pindel_dbsnp_filter.indel.input";
    print STDERR "Writing to $out\n";
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
pindel.dbsnp.indel.annotator = $snpsift_jar
pindel.dbsnp.indel.db = $dbsnp_db
pindel.dbsnp.indel.rawvcf = $filter_out
pindel.dbsnp.indel.mode = filter
pindel.dbsnp.indel.passfile  = $dbsnp_filtered_fn
pindel.dbsnp.indel.dbsnpfile = $filter_results/pindel.out.current_final.dbsnp_present.vcf
EOF

# 1. run pindel_filter.  This produces
#    pindel.out.raw.CvgVafStrand_pass 
#    pindel.out.raw.CvgVafStrand_fail
#    pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf  -> this is input into dbSnP filter
#    pindel.out.raw.CvgVafStrand_pass.Homopolymer_fail.vcf  
# 2. rename headers of pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf to be "NORMAL" and "TUMOR" 
#    Also add "FORMAT" column which new version of pindel2vcf misses
# 3. Run dbSnP filter
#    pindel.out.current_final.dbsnp_pass.vcf
#    pindel.out.current_final.dbsnp_pass.vcf.idx
#    pindel.out.current_final.dbsnp_present.vcf
# 4. run vcf_filter.py family of filters: VAF, read depth, and indel length
#    * Reads pindel.out.current_final.dbsnp_pass.vcf
#    * Outputs pindel.out.current_final.dbsnp_pass.filtered.vcf
# 5. Optionally delete intermediate files
#    - specifically, files with "_fail" in the filename

    my $vcf_filtered_fn = "$filter_results/pindel.out.current_final.dbsnp_pass.filtered.vcf";

    my $outfn = "$job_files_dir/$current_job_file";
    print STDERR "Writing to $outfn\n";
    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash

echo Running pindel_filter
$perl $gvip_dir/pindel_filter.pl $filter_results/pindel_filter.input
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

echo Running dbsnp_filter.pl
export JAVA_OPTS=\"-Xms256m -Xmx10g\"
$perl $gvip_dir/dbsnp_filter.pl $filter_results/pindel_dbsnp_filter.indel.input
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


echo Running combined vcf_filter.py filters: VAF, read depth, and indel length
export PYTHONPATH="$filter_dir:\$PYTHONPATH"
bash $filter_dir/run_combined_vcf_filter.sh $dbsnp_filtered_fn pindel $pindel_vcf_filter_config $vcf_filtered_fn
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


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
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print STDERR "Executing:\n $bsub_com \n";

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
