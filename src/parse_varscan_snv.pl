# Process varscan SNV output in 2 steps:
# * varscan processSomatic
# * varscan somaticFilter
#
# Note that dbSnP and length/depth/vaf filtering is removed from this step

# The following files created in $results_dir/varscan_out are read here:
#  varscan.out.som_indel.vcf 
#  varscan.out.som_snv.vcf 
#
# While the primary goal is to process the snv calls, indel calls are used as input for the somaticFilter
#
# processing which takes place here will be written to varscan/filter_snv_out 
#
# Principal output:
#   $results_dir/varscan/filter_snv_out/varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf   
#
# Note that all filenames above dependent on the filename of input data.

# The following parameters are read from varscan_config
#    snv.min-tumor-freq
#    snv.max-normal-freq
#    snv.p-value
#    filter.min-coverage
#    filter.min-reads2
#    filter.min-strands2
#    filter.min-avg-qual
#    filter.min-var-freq
#    filter.p-value

use File::Basename;

# Confirm that all required configuration parameters are defined.  Exit with an error if they are not
sub test_config_parameters_varscan_parse {
    my ($config_fn, %params) = @_;

    my @required_keys = (
        "snv.min-tumor-freq",
        "snv.max-normal-freq",
        "snv.p-value",
        "filter.min-coverage",
        "filter.min-reads2",
        "filter.min-strands2",
        "filter.min-avg-qual",
        "filter.min-var-freq",
        "filter.p-value");

    foreach my $key (@required_keys) {
        if (! exists $params{$key}) {
            die ("Required key $key not found in configuration file $config_fn\n");
        }
    }
}

# VarScan is unfortunate in that all output data is written to the same directory as input data, and
# the documentation does not describe a way to change that.  Since input data is passed, and we need
# to be able to control where data is written to, we must create a soft-link to input data in output 
# directory.  Note that link must be created with absolute, not relative, path
#
# make_data_link input_dat results_dir
#   creates a link to input_dat in results_dir
sub make_data_link {
    my $input_dat = shift;
    my $results_dir = shift;

    die "Error: Input file $input_dat does not exist\n" if (! -e $input_dat);
    $input_dat_abs = `readlink -f $input_dat`;
    chomp $input_dat_abs;
    print STDERR "Making link to $input_dat_abs from $results_dir\n";
    my $result = system ("ln -fs $input_dat_abs $results_dir "); 
    die("Exiting ($result).\n") if $result != 0;
    return ($results_dir . "/" . basename($input_dat_abs)) ;
}

sub parse_varscan_snv {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $varscan_indel_raw = shift; 
    my $varscan_snv_raw = shift;  
    my $varscan_config = shift;

    # define output directory: varscan/filter_snv_out
    my $filter_results_dir = "$results_dir/varscan/filter_snv_out";
    print STDERR "Filter results: $filter_results_dir\n";
    system("mkdir -p $filter_results_dir");

    # make input data available in working directory via soft links
    my $snv_raw = make_data_link($varscan_snv_raw, $filter_results_dir); 
    my $indel_raw = make_data_link($varscan_indel_raw, $filter_results_dir); 

    # Read configuration file into %params
    my %params = get_config_params($varscan_config, 1);
    test_config_parameters_varscan_parse($varscan_config, %params);

    #
    # ProcessSomatic SNV parameters and command 
    #
    my $somatic_snv_params="--min-tumor-freq $params{'snv.min-tumor-freq'} --max-normal-freq $params{'snv.max-normal-freq'} --p-value $params{'snv.p-value'}";  
    print STDERR "Somatic SNV Params:\n$somatic_snv_params\n";

    # processSomatic creates the following in the same directory as the input data
    #   varscan.out.som_snv.Somatic.hc.vcf      -> used for SNV SNP filter below 
    #   varscan.out.som_snv.Somatic.vcf        
    #   varscan.out.som_snv.LOH.hc.vcf         
    #   varscan.out.som_snv.LOH.vcf            
    #   varscan.out.som_snv.Germline.hc.vcf    
    #   varscan.out.som_snv.Germline.vcf       
    # all this based on assumption that snv_raw filename is varscan.out.som_snv.vcf
    # Filename processSomaticOut determined by `varscan processSomatic` based on $snv_raw
    my $processSomaticOut="$filter_results_dir/varscan.out.som_snv.Somatic.hc.vcf";  

    my $process_somatic_cmd = "java \${JAVA_OPTS} -jar $SWpaths::varscan_jar processSomatic $snv_raw $somatic_snv_params ";

    #
    # Somatic Filter parameters and command 
    # http://varscan.sourceforge.net/using-varscan.html#v2.3_somaticFilter
    #
    my $somatic_filter_params="--min-coverage $params{'filter.min-coverage'} --min-reads2 $params{'filter.min-reads2'} " .
        "--min-strands2 $params{'filter.min-strands2'} --min-avg-qual $params{'filter.min-avg-qual'} " . 
        "--min-var-freq $params{'filter.min-var-freq'} --p-value $params{'filter.p-value'}";
    print STDERR "Somatic Filter Params:\n$somatic_filter_params\n";
    # we define filename of somatic_filter_out :
    my $somatic_filter_out ="$filter_results_dir/varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf";
    my $somatic_filter_cmd = "java \${JAVA_OPTS} -jar $SWpaths::varscan_jar somaticFilter $processSomaticOut $somatic_filter_params --indel-file $indel_raw --output-file $somatic_filter_out ";

    #
    # Construct composite script
    #
    my $outfn = "$job_files_dir/j4a_parse_varscan_snv.sh";
    print STDERR "Writing to $outfn\n";
    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xms256m -Xmx10g\"

>&2 echo Applying varscan process filter to somatic SNVs
$process_somatic_cmd 
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

>&2 echo Applying varscan somatic filter
$somatic_filter_cmd
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

EOF

    close OUT;
    my $cmd = "bash < $outfn\n";
    print STDERR "Executing:\n $cmd \n";

    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
