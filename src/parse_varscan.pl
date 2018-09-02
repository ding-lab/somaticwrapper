# The following files were created in $sample_full_path/varscan_out
#  bamfilelist.inp
#  varscan.out.som.log
#  varscan.out.som_indel.vcf 
#  varscan.out.som_snv.vcf 
# processing which takes place here will be written to varscan/filter_out 
#
# Principal output and CWL mapping:
    # varscan.out.som_snv.Somatic.hc.vcf      -> varscan_snv_process
    # varscan.out.som_indel.Somatic.hc.vcf    -> varscan_indel_process
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf   -> varscan_snv_filtered
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.filtered.vcf   -> varscan_snv_dbsnp, used for merge_vcf
    # varscan.out.som_indel.Somatic.hc.dbsnp_pass.filtered.vcf    -> varscan_indel_dbsnp, used for merge_vcf

# The following parameters are read from varscan_config.  Numbers provided are parameters used by Song circa May 2018
#    snv.min-tumor-freq = 0.05 
#    snv.max-normal-freq = 0.05 
#    snv.p-value = 0.05
#    indel.min-tumor-freq = 0.05 
#    indel.max-normal-freq = 0.05 
#    indel.p-value = 0.05
#    filter.min-coverage = 20 
#    filter.min-reads2 = 4 
#    filter.min-strands2 = 1 
#    filter.min-avg-qual = 20 
#    filter.min-var-freq = 0.05 
#    filter.p-value = 0.05

# CWL changes:
# * get rid of genomevip_label steps
# * input filenames are passed explicitly: varscan_indel_raw and varscan_snv_raw
#   * these had been varscan.out.som_indel.vcf and varscan.out.som_snv.vcf

use File::Basename;

# Confirm that all required configuration parameters are defined.  Exit with an error if they are not
sub test_config_parameters_varscan_parse {
    my ($config_fn, %params) = @_;

    my @required_keys = (
        "snv.min-tumor-freq",
        "snv.max-normal-freq",
        "snv.p-value",
        "indel.min-tumor-freq",
        "indel.max-normal-freq",
        "indel.p-value",
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

sub parse_varscan{
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $filter_dir = shift;
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $varscan_jar = shift;
    my $varscan_indel_raw = shift; # indeloutgvip
    my $varscan_snv_raw = shift;  # snvoutgvip
    my $varscan_config = shift;
    my $varscan_vcf_filter_config = shift;

    $current_job_file = "j4_parse_varscan.sh";
    
    # It would be helpul to allow dbsnp_db to be not set, which would imply skipping the filtering step
    # This is currently not supported: require dbsnp_db to be defined and a file
    if ($dbsnp_db eq "") {
        die("Error: dbsnp_db not defined\n");
    } else {
        die "Error: dbSnP database file $dbsnp_db does not exist\n" if (! -e $dbsnp_db);
    }

    # Read configuration file into %params
    my %params = get_config_params($varscan_config, 1);
    test_config_parameters_varscan_parse($varscan_config, %params);

    # TODO: document these parameters
    # construct parameter arguments from %params
    my $somatic_snv_params="--min-tumor-freq $params{'snv.min-tumor-freq'} --max-normal-freq $params{'snv.max-normal-freq'} --p-value $params{'snv.p-value'}";  
    my $somatic_indel_params="--min-tumor-freq $params{'indel.min-tumor-freq'} " . 
        "--max-normal-freq $params{'indel.max-normal-freq'} --p-value $params{'indel.p-value'}";  
    my $somatic_filter_params="--min-coverage $params{'filter.min-coverage'} --min-reads2 $params{'filter.min-reads2'} " .
        "--min-strands2 $params{'filter.min-strands2'} --min-avg-qual $params{'filter.min-avg-qual'} " . 
        "--min-var-freq $params{'filter.min-var-freq'} --p-value $params{'filter.p-value'}";

    print "Somatic SNV Params:\n$somatic_snv_params\n";
    print "Somatic Indel Params:\n$somatic_indel_params\n";
    print "Somatic Filter Params:\n$somatic_filter_params\n";

    my $bsub = "bash";
    my $filter_results = "$sample_full_path/varscan/filter_out";
    print STDERR "Filter results: $filter_results\n";
    system("mkdir -p $filter_results");

    # VarScan is pathological in that all output data is written to the same directory as input data, and
    # the documentation does not describe a way to change that.  Since input data is passed, and we need
    # to be able to control where data is written to, we must create a soft-link to input data in output 
    # directory.  Note that link must be created with absolute, not relative, path
    die "Error: Indel raw input file $varscan_indel_raw does not exist\n" if (! -e $varscan_indel_raw);
    die "Error: SNV raw input file $varscan_snv_raw does not exist\n" if (! -e $varscan_snv_raw);

    $varscan_indel_raw = `readlink -f $varscan_indel_raw`;
    chomp $varscan_indel_raw;
    $varscan_snv_raw = `readlink -f $varscan_snv_raw`;
    chomp $varscan_snv_raw;

    system ("ln -fs $varscan_indel_raw $filter_results "); 
    system ("ln -fs $varscan_snv_raw $filter_results "); 
    my $indel_raw=$filter_results . "/" . basename($varscan_indel_raw) ;
    my $snv_raw=$filter_results . "/" . basename($varscan_snv_raw) ;

    # These based on original script

    my $thissnvorig="$filter_results/varscan.out.som_snv.Somatic.hc.vcf";  # This is genrated by varscan processSomatic
    my $myindelorig="$filter_results/varscan.out.som_indel.vcf";

    my $log_file="$filter_results/varscan.out.som.log";

    my $somsnvpass="$filter_results/varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf";

    # $dbsnp_filtered_snv_fn is the output of dbsnp filter of SNV calls
    my $dbsnp_filtered_snv_fn = "$filter_results/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf";
    my $out = "$filter_results/vs_dbsnp_filter.snv.input";
    print STDERR "Writing to $out\n";
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.snv.annotator = $snpsift_jar
varscan.dbsnp.snv.db = $dbsnp_db
varscan.dbsnp.snv.rawvcf = $somsnvpass
varscan.dbsnp.snv.mode = filter
varscan.dbsnp.snv.passfile  = $dbsnp_filtered_snv_fn
varscan.dbsnp.snv.dbsnpfile = $filter_results/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_present.vcf
EOF

    # $dbsnp_filtered_indel_fn is the output of dbsnp filter of SNV calls
    my $dbsnp_filtered_indel_fn = "$filter_results/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf";
    my $out = "$filter_results/vs_dbsnp_filter.indel.input";
    print STDERR "Writing to $out\n";
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.indel.annotator = $snpsift_jar
varscan.dbsnp.indel.db = $dbsnp_db
varscan.dbsnp.indel.rawvcf = $filter_results/varscan.out.som_indel.Somatic.hc.vcf
varscan.dbsnp.indel.mode = filter
varscan.dbsnp.indel.passfile  = $dbsnp_filtered_indel_fn
varscan.dbsnp.indel.dbsnpfile = $filter_results/varscan.out.som_indel.Somatic.hc.dbsnp_present.vcf
EOF

    # Run vcf_filter.py family of filters: VAF, read depth, and indel length
    #    * Reads varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf
    #        * Outputs varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.filtered.vcf
    #    * Reads varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf
    #        * Outputs varscan.out.som_indel.Somatic.hc.dbsnp_pass.filtered.vcf
    my $vcf_filtered_snv_fn = "$filter_results/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.filtered.vcf";
    my $vcf_filtered_indel_fn = "$filter_results/varscan.out.som_indel.Somatic.hc.dbsnp_pass.filtered.vcf";


    my $outfn = "$job_files_dir/$current_job_file";
    print STDERR "Writing to $outfn\n";
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xms256m -Xmx10g\"

echo \'APPLYING PROCESS FILTER TO SOMATIC SNVS:\' # &> $log_file
# Script below creates the following in the same directory as the input data
# The inability to define output directory complicates things
    # varscan.out.som_snv.Somatic.hc.vcf      -> used for SNV SNP filter below 
    # varscan.out.som_snv.Somatic.vcf        
    # varscan.out.som_snv.LOH.hc.vcf         
    # varscan.out.som_snv.LOH.vcf            
    # varscan.out.som_snv.Germline.hc.vcf    
    # varscan.out.som_snv.Germline.vcf       
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $snv_raw $somatic_snv_params /varscan # &>> $log_file
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


echo \'APPLYING PROCESS FILTER TO SOMATIC INDELS:\' # &>> $log_file
# Script below creates:
    # varscan.out.som_indel.Germline.hc.vcf    
    # varscan.out.som_indel.Germline.vcf       
    # varscan.out.som_indel.LOH.hc.vcf         
    # varscan.out.som_indel.LOH.vcf            
    # varscan.out.som_indel.Somatic.hc.vcf     -> used for Indel SnP Filter below 
    # varscan.out.som_indel.Somatic.vcf        
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $indel_raw   $somatic_indel_params  # &>> $log_file
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi



### Somatic Filter filters SNV based on indel
# http://varscan.sourceforge.net/using-varscan.html#v2.3_somaticFilter
echo \'APPLYING SOMATIC FILTER:\' # &>> $log_file

# Script below creates:
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf   -> used for SNV dbSnP 
java \${JAVA_OPTS} -jar $varscan_jar somaticFilter  $thissnvorig $somatic_filter_params  --indel-file  $indel_raw --output-file  $somsnvpass  # &>> $log_file   
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


### dbSnP Filter

# 1) SNV
# Script below reads:  
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf
# and generates:
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_present.vcf  
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf     -> used for vcf_filter.py
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_anno.vcf   
$perl $gvip_dir/dbsnp_filter.pl  $filter_results/vs_dbsnp_filter.snv.input
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


# 2) indel
# Script below reads
    # varscan.out.som_indel.Somatic.hc.vcf
# and generates:
    # varscan.out.som_indel.Somatic.hc.dbsnp_present.vcf 
    # varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf    -> used for vcf_filter.py
    # varscan.out.som_indel.Somatic.hc.dbsnp_anno.vcf    
$perl $gvip_dir/dbsnp_filter.pl $filter_results/vs_dbsnp_filter.indel.input
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi


    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.filtered.vcf     -> used for merge_vcf
    # varscan.out.som_indel.Somatic.hc.dbsnp_pass.filtered.vcf    -> used for merge_vcf

echo Running combined vcf_filter.py filters: VAF, read depth, and indel length
export PYTHONPATH="$filter_dir:\$PYTHONPATH"
bash $filter_dir/run_combined_vcf_filter.sh $dbsnp_filtered_snv_fn varscan $varscan_vcf_filter_config $vcf_filtered_snv_fn 
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi

bash $filter_dir/run_combined_vcf_filter.sh $dbsnp_filtered_indel_fn varscan $varscan_vcf_filter_config $vcf_filtered_indel_fn 
rc=\$?
if [[ \$rc != 0 ]]; then
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc;
fi



EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print STDERR "Executing:\n $bsub_com \n";

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
