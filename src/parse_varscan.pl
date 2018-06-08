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
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf   -> varscan_snv_dbsnp
    # varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf    -> varscan_indel_dbsnp

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
    my %params = shift;
    my $config_fn = shift;

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
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $varscan_jar = shift;
    my $varscan_indel_raw = shift; # indeloutgvip
    my $varscan_snv_raw = shift;  # snvoutgvip
    my $varscan_config = shift;

    $current_job_file = "j4_parse_varscan.sh";
    die "Error: dbSnP database file $dbsnp_db does not exist\n" if (! -e $dbsnp_db);

    # Read configuration file into %params
    my %params = get_config_params($varscan_config, 1);
    test_config_parameters_varscan_parse(%params, $varscan_config);

    # TODO: document these parameters
    # construct parameter arguments from %params
    my $somatic_snv_params="--min-tumor-freq $params{'snv.min-tumor-freq'} --max-normal-freq $params{'snv.max-normal-freq'} --p-value $params{'snv.p-value'}";  
    my $somatic_indel_params="--min-tumor-freq $params{'indel.min-tumor-freq'} --max-normal-freq $params{'indel.max-normal-freq'} --p-value $params{'indel.p-value'}";  
    my $somatic_filter_params="--min-coverage $params{'filter.min-coverage'} --min-reads2 $params{'filter.min-reads2'} " +
        "--min-strands2 $params{'filter.min-strands2'} --min-avg-qual $params{'filter.min-avg-qual'} " + 
        "--min-var-freq $params{'filter.min-var-freq'} --p-value $params{'filter.p-value'}";

    print "Somatic SNV Params:\n$somatic_snv_params\n";
    print "Somatic Indel Params:\n$somatic_indel_params\n";
    print "Somatic Filter Params:\n$somatic_filter_params\n";
die("Quitting early\n");

    my $bsub = "bash";
    my $filter_results = "$sample_full_path/varscan/filter_out";
    system("mkdir -p $filter_results");

    # VarScan is pathological in that all output data is written to the same directory as input data, and
    # the documentation does not describe a way to change that.  Since input data is passed, and we need
    # to be able to control where data is written to, we must create a soft-link to input data in output 
    # directory.
    die "Error: Indel raw input file $varscan_indel_raw does not exist\n" if (! -e $varscan_indel_raw);
    die "Error: SNV raw input file $varscan_snv_raw does not exist\n" if (! -e $varscan_snv_raw);

    system ("ln -fs $varscan_indel_raw $filter_results "); 
    system ("ln -fs $varscan_snv_raw $filter_results "); 
    my $indel_raw=$filter_results . "/" . basename($varscan_indel_raw) ;
    my $snv_raw=$filter_results . "/" . basename($varscan_snv_raw) ;

    # These based on original script
    my $snvoutbase="$filter_results/varscan.out.som_snv";
    my $indeloutbase="$filter_results/varscan.out.som_indel";

    #my $thissnvorig="${snvoutbase}.Somatic.hc.vcf";  # This is genrated by varscan processSomatic
    my $thissnvorig="$filter_results/varscan.out.som_snv.Somatic.hc.vcf";  # This is genrated by varscan processSomatic
    my $myindelorig="${indeloutbase}.vcf";
    my $thissnvpass="${snvoutbase}.Somatic.hc.somfilter_pass.vcf";


    my $log_file="$filter_results/varscan.out.som.log";

    my $somsnvpass="$filter_results/varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf";

    my $indeloutbase="$filter_results/varscan.out.som_indel";

    my $out = "$filter_results/vs_dbsnp_filter.snv.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.snv.annotator = $snpsift_jar
varscan.dbsnp.snv.db = $dbsnp_db
varscan.dbsnp.snv.rawvcf = $somsnvpass
varscan.dbsnp.snv.mode = filter
varscan.dbsnp.snv.passfile  = $filter_results/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf
varscan.dbsnp.snv.dbsnpfile = $filter_results/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_present.vcf
EOF

    my $out = "$filter_results/vs_dbsnp_filter.indel.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.indel.annotator = $snpsift_jar
varscan.dbsnp.indel.db = $dbsnp_db
varscan.dbsnp.indel.rawvcf = $indeloutbase.Somatic.hc.vcf
varscan.dbsnp.indel.mode = filter
varscan.dbsnp.indel.passfile  = $indeloutbase.Somatic.hc.dbsnp_pass.vcf
varscan.dbsnp.indel.dbsnpfile = $indeloutbase.Somatic.hc.dbsnp_present.vcf
EOF

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xms256m -Xmx10g\"

echo \'APPLYING PROCESS FILTER TO SOMATIC SNVS:\' # &> $log_file
# Script below creates the following in the same directory as the input data
# The inability to define output directory is maddening 
    # varscan.out.som_snv.Somatic.hc.vcf      -> used for SNV SNP filter below and vep annotation
    # varscan.out.som_snv.Somatic.vcf        
    # varscan.out.som_snv.LOH.hc.vcf         
    # varscan.out.som_snv.LOH.vcf            
    # varscan.out.som_snv.Germline.hc.vcf    
    # varscan.out.som_snv.Germline.vcf       
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $snv_raw $somatic_snv_params /varscan # &>> $log_file

echo \'APPLYING PROCESS FILTER TO SOMATIC INDELS:\' # &>> $log_file
# Script below creates:
    # varscan.out.som_indel.Germline.hc.vcf    
    # varscan.out.som_indel.Germline.vcf       
    # varscan.out.som_indel.LOH.hc.vcf         
    # varscan.out.som_indel.LOH.vcf            
    # varscan.out.som_indel.Somatic.hc.vcf     -> used for Indel SnP Filter below and vep annotation
    # varscan.out.som_indel.Somatic.vcf        
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $indel_raw   $somatic_indel_params  # &>> $log_file


### Somatic Filter filters SNV based on indel
# http://varscan.sourceforge.net/using-varscan.html#v2.3_somaticFilter
echo \'APPLYING SOMATIC FILTER:\' # &>> $log_file

# Script below creates:
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf   -> used for SNV dbSnP and vep annotation
java \${JAVA_OPTS} -jar $varscan_jar somaticFilter  $thissnvorig $somatic_filter_params  --indel-file  $indel_raw --output-file  $somsnvpass  # &>> $log_file   

### dbSnP Filter

# 1) SNV
# Script below reads:  
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.vcf
# and generates:
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_present.vcf  
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf     -> used for merge_vcf
    # varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_anno.vcf   
$perl $gvip_dir/dbsnp_filter.pl  $filter_results/vs_dbsnp_filter.snv.input

# 2) indel
# Script below reads
    # varscan.out.som_indel.Somatic.hc.vcf
# and generates:
    # varscan.out.som_indel.Somatic.hc.dbsnp_present.vcf 
    # varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf    -> used for merge_vcf
    # varscan.out.som_indel.Somatic.hc.dbsnp_anno.vcf    
$perl $gvip_dir/dbsnp_filter.pl $filter_results/vs_dbsnp_filter.indel.input

# Genome VIP SNV filter , and two GenomeVIP VEP annotation calls deleted from end of workflow per 
# discussion with Song

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

1;
