# Annotate VCF is a CWL successor to run_vep.  It converts a given VCF to VEP
# format (or just annotates vcf).  It operates on just one input file, however,
# and writes to a given output file.

# if $cache_dir is defined, use that as path to VEP Cache
# if $cache_gz is defined, then copy that .tar.gz file to $cache_dir and extract it there
#   (this is used for a cwl setup where arbitrary paths are not accessible)

# Logic of use_vep_db is defined by $cache_dir.  If $cache_dir is set to a value, then use_vep_db = 0.
#   * also test to see if is a directory.  
# if use_vep_db is 1, run step in db mode.
#   this 1) uses online database (so cache isn't installed) 2) does not use tmp files
#   It is meant to be used for testing and lightweight applications.  Use the cache for
#   better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 

# assembly is the assembly argument passed to vep
# output_vep is a boolean.  Output annotated VEP rather than VCF format.
# cache_dir defined implies use_vep_db is 0
#   Cache installation is done in somaticwrapper/image.setup/D_VEP
#   See https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
# 
# if $cache_gz is defined, it is assumed this is a .tar.gz version of VEP cache.
#   extract its contents into $cache_dir (./vep-cache if not specified)
#   It will subsequently be deleted

# helper function
sub write_vep_input {
    my $config_fn = shift;
    my $module = shift;  # e.g. varscan.vep or strelka.vep
    my $vcf = shift;
    my $output = shift; 
    my $vep_cmd = shift;
    my $cache_dir = shift;   
    my $REF = shift;
    my $assembly = shift;
    my $cache_version = shift;   
    my $use_vep_db = shift;  # 1 for testing/demo, 0 for production
    my $output_vep = shift;  # output annotated vep rather than vcf format after merge step.  add suffix 'vep' to output

    if ($output_vep) {
        $output = "$output.vep";
    }

    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$module.vcf = $vcf
$module.output = $output
$module.vep_cmd = $vep_cmd
$module.cachedir = $cache_dir
$module.reffasta = $REF
$module.assembly = $assembly
$module.cache_version = $cache_version
$module.usedb = $use_vep_db  
$module.output_vep = $output_vep  
EOF
}

sub annotate_vcf {
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $assembly = shift;
    my $cache_version = shift; # 90
    my $cache_dir = shift;  # if defined, implies use_vep_db = 0
    my $cache_gz = shift;   
    my $output_vep = shift;  # if 1, output annotated vep after merge step.  If 0, output vcf format 
    my $input_vcf = shift;  # for CWL work, we are passed an input VCF
    my $output_vcf = shift;  # for CWL work, we are passed an output annotated vcf

    $current_job_file = "j10_vep.sh";

    my $bsub = "bash";
    my $filter_results = "$sample_full_path/vep";
    system("mkdir -p $filter_results");

    my $config_fn = "$filter_results/vep.merged.input";

    my $use_vep_db = 1;

    # if cache_gz is defined, extract its contents to $cache_dir
    # if $cache_dir is not defined, defaults to ./vep-cache
    if ( $cache_gz ) {
        if (! $cache_dir) {
            $cache_dir = "./vep-cache";
        }
        if (! -d $cache_dir) {
            mkdir $cache_dir or die "$!\n";
        }
        print STDERR "Extracting VEP Cache tarball $cache_gz into $cache_dir";
        my $rc = system ("tar -zxf $cache_gz --directory $cache_dir");
        die("Exiting ($rc).\n") if $rc != 0;
    }

    if ( $cache_dir ) {
        die "\nError: Cache dir $cache_dir does not exist\n" if (! -d $cache_dir);
        die "\nError: Please specify --cache_version \n" if (! $cache_version);
        $use_vep_db = 0;
    }

# Adding final output to vep annotation
    write_vep_input(
        $config_fn,          # Config fn
        "merged.vep",                               # Module
        $input_vcf,                # VCF (input)
        $output_vcf,           # output
        $vep_cmd, $cache_dir, $REF, $assembly, $cache_version, $use_vep_db, $output_vep);

    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx512m\"

$perl $gvip_dir/vep_annotator.pl $config_fn

# Evaluate return value see https://stackoverflow.com/questions/90418/exit-shell-script-based-on-process-exit-code
rc=\$? 
if [[ \$rc != 0 ]]; then 
    >&2 echo Fatal error \$rc: \$!.  Exiting.
    exit \$rc; 
fi

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;

    # Clean up by deleting contents of cache_dir - this tends to be big (>10Gb) and unnecessary to keep
    if ( $cache_gz ) {
        print STDERR "Deleting $cache_dir\n";
        my $rc = system("rm -rf $cache_dir\n");
        die("Exiting ($rc).\n") if $rc != 0;
    }
}

1;
