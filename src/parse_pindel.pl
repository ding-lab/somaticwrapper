
# Skipping VEP annotation

# principal output used in merging: pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf
# -> pindel_dbsnp output

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
    my $pindel_dir = shift;
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $pindel_config = shift;
    my $pindel_raw_in = shift; # NEW
    my $no_delete_temp = shift;

    if (! $no_delete_temp) {
        $no_delete_temp = 0; # avoid empty variables
    }

    $current_job_file = "j7_parse_pindel.sh";

    my $bsub = "bash";
    my $filter_results = "$sample_full_path/pindel/filter_out";
    system("mkdir -p $filter_results");

    die "Error: dbSnP database file $dbsnp_db does not exist\n" if (! -e $dbsnp_db);

    # pindel_filter is pathological in that all output data is written to the same directory as input data, and
    # the documentation does not describe a way to change that.  Since input data is passed, and we need
    # to be able to control where data is written to, we must create a soft-link to input data in output 
    # directory.  Link has to be created with an absolute filename
    die "Error: Pindel raw input file $pindel_raw_in does not exist\n" if (! -e $pindel_raw_in);
    $pindel_raw_in = `readlink -f $pindel_raw_in`;
    chomp $pindel_raw_in;

    system ("ln -fs $pindel_raw_in $filter_results "); 
    my $pindel_raw=$filter_results . "/" . basename($pindel_raw_in) ;

    my $outlist="$filter_results/pindel.out.filelist";
    #my $current_final="$filter_results/pindel.out.current_final.Somatic.vcf";
    # This is the principal result of pindel_filter
    my $filter_out="$pindel_raw.CvgVafStrand_pass.Homopolymer_pass.vcf";

## Pindel Filter - below is input into pindel_filter.v0.5
# lines below are added to data from $pindel_config
    die "$pindel_config does not exist\n" unless (-f $pindel_config);

    my $out = "$filter_results/pindel_filter.input";
    print("Copying $pindel_config to $out and appending\n");
    system("cp $pindel_config $out");

    open(OUT, ">>$out") or die $!;
    print OUT <<"EOF";
pindel.filter.pindel2vcf = $pindel_dir/pindel2vcf
pindel.filter.variants_file = $pindel_raw
pindel.filter.REF = $REF
pindel.filter.date = 000000
EOF

## dbSnP Filter
    my $out = "$filter_results/pindel_dbsnp_filter.indel.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
pindel.dbsnp.indel.annotator = $snpsift_jar
pindel.dbsnp.indel.db = $dbsnp_db
pindel.dbsnp.indel.rawvcf = $filter_out
pindel.dbsnp.indel.mode = filter
pindel.dbsnp.indel.passfile  = $filter_results/pindel.out.current_final.dbsnp_pass.vcf
pindel.dbsnp.indel.dbsnpfile = $filter_results/pindel.out.current_final.dbsnp_present.vcf
EOF

# 1. run pindel_filter.  This produces
#    pindel.out.raw.CvgVafStrand_pass 
#    pindel.out.raw.CvgVafStrand_fail
#    pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf  -> this is input into dbSnP filter
#    pindel.out.raw.CvgVafStrand_pass.Homopolymer_fail.vcf  
# 2. rename headers of pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf to be "NORMAL" and "TUMOR" 
# 3. Run dbSnP filter
#    pindel.out.current_final.dbsnp_pass.vcf
#    pindel.out.current_final.dbsnp_pass.vcf.idx
#    pindel.out.current_final.dbsnp_present.vcf
# 4. Optionally delete intermediate files
#    - specifically, files with "_fail" in the filename

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash

echo Running pindel_filter.v0.5.pl
$perl $gvip_dir/pindel_filter.v0.5.pl $filter_results/pindel_filter.input

# Reheader output of pindel_filter to have sample names "NORMAL" and "TUMOR"
TMP=$filter_out.tmp
mv $filter_out $TMP
awk 'BEGIN{FS="\t";OFS="\t"}{if ($1 == "#CHROM") print $1, $2, $3, $4, $5, $6, $7, $8, $9, "NORMAL", "TUMOR"; else print}' $TMP > $filter_out

export JAVA_OPTS=\"-Xms256m -Xmx10g\"

echo Running dbsnp_filter.pl
$perl $gvip_dir/dbsnp_filter.pl $filter_results/pindel_dbsnp_filter.indel.input

if [[ $no_delete_temp == 1 ]]; then

>&2 echo Not deleting intermediate files

else

>&2 echo Deleting intermediate \\"filter fail\\" files
cd $filter_results
rm -f \*_fail\* $TMP

fi

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
