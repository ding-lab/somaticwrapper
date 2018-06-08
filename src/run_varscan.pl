# Generates the files in varscan/varscan_out
# * bamfilelist.inp
# * varscan.out.som_indel.vcf
# * varscan.out.som_snv.vcf

# return hash of parameters of form params{key}=value, where key and value are specified in configuration file as
#   key = value
# Unlike params in GenomeVIP, where a key of the form "x.y.z" is stripped so key=z, here the entire key is retained
sub get_config_params {
    my $config_fn = shift;
    my $DEBUG=shift;

    print ("Reading configuration file: $config_fn\n");

    open( my $fh, '<', $config_fn ) or die "Can't open config file $config_fn: $!";

    my %paras;
    # first form is from GenomeVIP/dbsnp_filter.pl
    # map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
    map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ $_[0] } = $v } } (<$fh>);
    close $fh;

    if ($DEBUG) {
        map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;
    }
    return %paras;
}

# Confirm that all required configuration parameters are defined.  Exit with an error if they are not
sub test_config_parameters_varscan_run {

# Discussion of subtleties of passing hashes in perl: pass strings first
# https://stackoverflow.com/questions/3567316/unable-to-pass-a-hash-and-a-string-to-a-function-together-in-perl
    
    my ($config_fn, %params) = @_;

    my @required_keys = (
        "varscan.mpileup",
        "varscan.p-value",
        "varscan.somatic-p-value",
        "varscan.min-coverage-normal",
        "varscan.min-coverage-tumor",
        "varscan.min-var-freq",
        "varscan.min-freq-for-hom",
        "varscan.normal-purity",
        "varscan.tumor-purity",
        "varscan.strand-filter",
        "varscan.min-avg-qual");

    foreach my $key (@required_keys) {
        if (! exists $params{$key}) {
            die ("Required key $key not found in configuration file $config_fn\n");
        }
    }
}

sub run_varscan{
    my $IN_bam_T = shift;
    my $IN_bam_N = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $varscan_config = shift;
    my $varscan_jar = shift;

    my $bsub = "bash";
    my $samtools="/usr/local/bin/samtools";

    $current_job_file = "j2_varscan.sh";
    die "Error: Tumor BAM $IN_bam_T does not exist\n" if (! -e $IN_bam_T);
    die "Error: Tumor BAM $IN_bam_T is empty\n" if (! -s $IN_bam_T);
    die "Error: Normal BAM $IN_bam_N does not exist\n" if (! -e $IN_bam_N);
    die "Error: Normal BAM $IN_bam_N is empty\n" if (! -s $IN_bam_N);

    my $workdir="$sample_full_path/varscan/varscan_out";
    system("mkdir -p $workdir");

    # Create a list of BAM files for varscan to use
    my $bam_list="$workdir/bamfilelist.inp";
    open(OUT, ">$bam_list") or die $!;
    print OUT "$IN_bam_N\n";
    print OUT "$IN_bam_T\n"; 
    close OUT;


    # Create the run script
    # Using HERE docs: https://stackoverflow.com/questions/17479354/how-to-use-here-doc-to-print-lines-to-file

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;


    my $run_name="varscan.out.som";
    my $log=$workdir."/".$run_name.".log";
    my $snvout=$workdir."/".$run_name."_snv";
    my $indelout=$workdir."/".$run_name."_indel";

#    die "File not found: $varscan_config\n" if (! -e $varscan_config);
#    # ignore comments in varscan_config and convert newlines to spaces, so that all arguments are in one line
#    my $varscan_args=`grep -v "^#" $varscan_config | tr '\n' ' '`;
#
#--mpileup 1 
#--p-value 0.99 
#--somatic-p-value 0.05 
#--min-coverage-normal 20 
#--min-coverage-tumor 20 
#--min-var-freq 0.05 
#--min-freq-for-hom 0.75 
#--normal-purity 1.00 
#--tumor-purity 1.00 
#--strand-filter 1 
#--min-avg-qual 15 

    # Read configuration file into %params
    my %params = get_config_params($varscan_config, 0);

    test_config_parameters_varscan_run($varscan_config, %params);

   my $varscan_args = "
--mpileup $params{'varscan.mpileup'} 
--p-value $params{'varscan.p-value'} 
--somatic-p-value $params{'varscan.somatic-p-value'} 
--min-coverage-normal $params{'varscan.min-coverage-normal'} 
--min-coverage-tumor $params{'varscan.min-coverage-tumor'} 
--min-var-freq $params{'varscan.min-var-freq'} 
--min-freq-for-hom $params{'varscan.min-freq-for-hom'} 
--normal-purity $params{'varscan.normal-purity'} 
--tumor-purity $params{'varscan.tumor-purity'} 
--strand-filter $params{'varscan.strand-filter'} 
--min-avg-qual $params{'varscan.min-avg-qual'} 
--output-vcf 1";

    $varscan_args =~ tr/\r\n//d;

    print OUT <<"EOF";
#!/bin/bash
JAVA_OPTS="-Xms256m -Xmx512m"

#echo Log to $log

SAMTOOLS_CMD="$samtools mpileup -q 1 -Q 13 -B -f $REF -b $bam_list "

JAVA_CMD="java \$JAVA_OPTS -jar $varscan_jar somatic - $run_name $varscan_args --output-snp $snvout --output-indel $indelout"

\$SAMTOOLS_CMD | \$JAVA_CMD # &> $log

EOF
    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";

    print("Executing:\n $bsub_com \n");

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;

}

1;
