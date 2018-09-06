# Process varscan indel output in 2 steps:
# * varscan processSomatic 
# * vcf_filter: VAF, Length, Depth
#
# Note that dbSnP filtering is removed from this step

# The following file created in $sample_full_path/varscan_out is read here:
#  varscan.out.som_indel.vcf 
#
# processing which takes place here will be written to varscan/filter_indel_out 
#
# Principal output 
    # varscan.out.som_indel.Somatic.hc.vcf    -> varscan_indel_process
    # varscan.out.som_indel.Somatic.hc.filtered.vcf    -> varscan_indel_dbsnp, used for merge_vcf

# The following parameters are read from varscan_config.  Numbers provided are parameters used by Song circa May 2018
#    indel.min-tumor-freq = 0.05 
#    indel.max-normal-freq = 0.05 
#    indel.p-value = 0.05

# * While input filenames are passed explicitly, internal naming logic still assumes that input data are named 
#   varscan.out.som_indel.vcf and varscan.out.som_snv.vcf

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
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $filter_dir = shift;
    my $varscan_jar = shift;
    my $varscan_indel_raw = shift; # indeloutgvip
    my $varscan_config = shift;
    my $varscan_vcf_filter_config = shift;
    
    my $filter_results = "$sample_full_path/varscan/filter_indel_out";
    print STDERR "Filter results: $filter_results\n";
    system("mkdir -p $filter_results");

    # creating a link to indel_raw.  See parse_varscan_snv.pl for details
    my $indel_raw = make_data_link($varscan_indel_raw, $filter_results_dir); 

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
        # varscan.out.som_indel.Somatic.hc.vcf     -> used for Indel SnP Filter below 
        # varscan.out.som_indel.Somatic.vcf        
    # all this based on assumption that indel_raw filename is varscan.out.som_indel.vcf
    my $process_somatic_out ="$filter_results_dir/varscan.out.som_indel.Somatic.hc.vcf";  
    my $process_somatic_cmd = "java \${JAVA_OPTS} -jar $varscan_jar processSomatic $indel_raw $somatic_indel_params";

    #
    # vcf_filter.py family of filters: VAF, read depth, and indel length
    #
    # we define filename of vcfFilteredSNVOut:
    my $vcf_filter_out = "$filter_results/varscan.out.som_indel.Somatic.hc.filtered.vcf";
    my $vcf_filter_cmd="bash $filter_dir/run_combined_vcf_filter.sh $process_somatic_out varscan $varscan_vcf_filter_config $vcf_filter_out ";

    my $current_job_file = "j4b_parse_varscan_indel.sh";
    my $outfn = "$job_files_dir/$current_job_file";
    print STDERR "Writing to $outfn\n";
    open(OUT, ">$outfn") or die $!;

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

>&2 echo Running combined vcf_filter.py filters: VAF, read depth, and indel length
export PYTHONPATH="$filter_dir:\$PYTHONPATH"
$vcf_filter_cmd
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
