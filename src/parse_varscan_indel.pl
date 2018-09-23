# Process varscan indel output using
# * varscan processSomatic 
#
# Note that dbSnP and length/depth/vaf filtering is removed from this step

# The following file created in $results_dir/varscan_out is read here:
#  varscan.out.som_indel.vcf 
#
# processing which takes place here will be written to varscan/filter_indel_out 
#
# Principal output 
#   $results_dir/varscan/filter_indel_out/varscan.out.som_indel.Somatic.hc.vcf    

# The following parameters are read from varscan_config.  
#    indel.min-tumor-freq
#    indel.max-normal-freq
#    indel.p-value

# * While input filenames are passed explicitly, internal naming logic still assumes that input data are named 
#   varscan.out.som_indel.vcf 

use File::Basename;

# Confirm that all required configuration parameters are defined.  Exit with an error if they are not
sub test_config_parameters_varscan_parse {
    my ($config_fn, %params) = @_;

    my @required_keys = (
        "indel.min-tumor-freq",
        "indel.max-normal-freq",
        "indel.p-value");

    foreach my $key (@required_keys) {
        if (! exists $params{$key}) {
            die ("Required key $key not found in configuration file $config_fn\n");
        }
    }
}

sub parse_varscan_indel {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $varscan_indel_raw = shift; 
    my $varscan_config = shift;
    
    my $filter_results = "$results_dir/varscan/filter_indel_out";
    print STDERR "Filter results: $filter_results\n";
    system("mkdir -p $filter_results");

    # creating a link to indel_raw.  See parse_varscan_snv.pl for details
    my $indel_raw = make_data_link($varscan_indel_raw, $filter_results); 

    # Read configuration file into %params
    my %params = get_config_params($varscan_config, 1);
    test_config_parameters_varscan_parse($varscan_config, %params);

    #
    # ProcessSomatic indel parameters and command 
    #
    # construct parameter arguments from %params
    my $somatic_indel_params="--min-tumor-freq $params{'indel.min-tumor-freq'} " . 
        "--max-normal-freq $params{'indel.max-normal-freq'} --p-value $params{'indel.p-value'}";  
    print STDERR "Somatic Indel Params:\n$somatic_indel_params\n";

    # processSomatic creates:
        # varscan.out.som_indel.Germline.hc.vcf    
        # varscan.out.som_indel.Germline.vcf       
        # varscan.out.som_indel.LOH.hc.vcf         
        # varscan.out.som_indel.LOH.vcf            
        # varscan.out.som_indel.Somatic.hc.vcf     -> principal output
        # varscan.out.som_indel.Somatic.vcf        
    # all this based on assumption that indel_raw filename is varscan.out.som_indel.vcf
    my $process_somatic_out ="$filter_results/varscan.out.som_indel.Somatic.hc.vcf";  
    my $process_somatic_cmd = "java \${JAVA_OPTS} -jar $SWpaths::varscan_jar processSomatic $indel_raw $somatic_indel_params";

    my $runfn = "$job_files_dir/j4b_parse_varscan_indel.sh";
    print STDERR "Writing to $runfn\n";
    open(OUT, ">$runfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xms256m -Xmx10g\"

>&2 echo Applying process filter to somatic indels
$process_somatic_cmd
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
