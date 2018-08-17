# Annotate VCF file, and write output as VCF, VEP, or MAF

# Annotation is performed using vep directly (if vep_output is 'vcf' or 'vep'), or vcf2maf (vep_output is 'maf'); the
# latter uses vep as well.  We can use either a local cache or online ('db') lookups (the latter implemented only for vep_output='vcf' or 'vep')..

#    --vep_cache_dir s: defines location of VEP cache directory and indicates whether to use online VEP DB lookups.  
#        * if vep_cache_dir is not defined, will perform online VEP DB lookups
#        * If vep_cache_dir is a directory, it indicates location of VEP cache 
#        * If vep_cache_dir is a file ending in .tar.gz, will extract its contents into "./vep-cache" and use VEP cache
#        NOTE: Online VEP database lookups a) uses online database (so cache isn't installed) b) does not use tmp files
#          It is meant to be used for testing and lightweight applications.  Use the cache for better performance.
#          See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
#
# if $vep_cache_dir is a .tar.gz file, then copy it to $cache_dir="./vep-cache" and extract it there
#   (assuming it is a .tar.gz version of VEP cache; this is typically used for a cwl setup where arbitrary paths are not accessible)
#   These contents will subsequently be deleted
#   Cache installation is done in somaticwrapper/image.setup/D_VEP
#   See https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
#
# assembly is the assembly argument passed to vep.  Optional
# vep_output: Defines output format after annotation.  Allowed values: vcf, vep, maf.  Default is vcf
# 
# Output is $results_dir/vep/output.vcf

my $perl = "/usr/bin/perl";  # this is for docker environment

sub annotate_vcf {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $reference = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $assembly = shift;
    my $cache_version = shift; # e.g., 90
    my $cache_dir = shift;  
    my $vep_output = shift;  # Output format following vep annotation.  May be 'vcf', 'vep', 'maf'
    my $input_vcf = shift;  # Name of input VCF to process

    # assembly and cache_version may be blank; if so, not passed on command line to vep

    if ($vep_output) {
        my @known_formats = ('vcf', 'vep', 'maf');
        my $format_OK = 0;
        foreach my $format (@known_formats) {
            if ($vep_output =~ $format) { $format_OK = 1; break }
        }
        if (not $format_OK) { die("Error: unknown VEP output format (--vep_output): $vep_output\n");}
    } else {
        $vep_output = "vcf";
    }


    $current_job_file = "j10_vep.sh";

    my $bsub = "bash";
    my $filter_results = "$results_dir/vep";
    system("mkdir -p $filter_results");

    my $config_fn = "$filter_results/vep.merged.input";
    my $output_fn;

    my $use_vep_db = 1;
    my $cache_gz;

    # if $cache_dir is a .tar.gz file, extract its contents to ./vep-cache
    # if $cache_dir is defined, confirm it exists
    # Otherwise, use VEP DB mode
    if ( $cache_dir =~ /\.tar\.gz/ ) {
        $cache_gz = $cache_dir;
        $cache_dir = "./vep-cache";
        if (! -d $cache_dir) {
            mkdir $cache_dir or die "$!\n";
        }
        print STDERR "Extracting VEP Cache tarball $cache_gz into $cache_dir\n";
        my $rc = system ("tar -zxf $cache_gz --directory $cache_dir");
        die("Exiting ($rc).\n") if $rc != 0;
        $use_vep_db = 0;
    } elsif (defined $cache_dir) {
        die "\nError: Cache dir $cache_dir does not exist\n" if (! -d $cache_dir);
        $use_vep_db = 0;
    }    

    # now need to decide if using vcf2maf or vep_annotator.pl
    my $cmd;
    if ($vep_output =~ /maf/) {
        $output_fn = "$filter_results/output.maf";

#      --ncbi-build     NCBI reference assembly of variants MAF (e.g. GRCm38 for mouse) [GRCh37]
#      --cache-version  Version of offline cache to use with VEP (e.g. 75, 84, 91) [Default: Installed version]
#      --vep-data       VEP's base cache/plugin directory [~/.vep]
        die ("--cache_dir must be defined for \"--vep_output maf\"\n") if !($cache_dir);

        my $opts = "--vep-data $cache_dir";
        if ($assembly) { $opts = "$opts --ncbi-build $assembly"; }
        if ($cache_version)  { $opts = "$opts --cache-version $cache_version"; }

        my $vep_path = dirname($vep_cmd);
        $cmd = "$perl /usr/local/mskcc-vcf2maf/vcf2maf.pl $opts --input-vcf $input_vcf --output-maf $output_fn --ref-fasta $reference --filter-vcf 0 --vep-path $vep_path --tmp-dir $filter_results";

    } else {
        if ($vep_output =~ /vcf/) {
            $output_fn = "$filter_results/output.vcf";
        } else {
            $output_fn = "$filter_results/output.vep";
        }

        # Adding final output to vep annotation
        write_vep_input(
            $config_fn,          # Config fn
            "merged.vep",                               # Module
            $input_vcf,                # VCF (input)
            $output_fn,           # output
            $vep_cmd, $cache_dir, $reference, $assembly, $cache_version, $use_vep_db, $vep_output =~ /vep/);

        $cmd = <<"EOF" 
        export JAVA_OPTS=\"-Xms256m -Xmx512m\"
        $perl $gvip_dir/vep_annotator.pl $config_fn
EOF
    }

    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

$cmd

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
    print STDERR "Final results written to $output_fn\n";
}

# helper function
sub write_vep_input {
    my $config_fn = shift;
    my $module = shift;  # e.g. varscan.vep or strelka.vep
    my $vcf = shift;
    my $output_fn = shift; 
    my $vep_cmd = shift;
    my $cache_dir = shift;   
    my $reference = shift;
    my $assembly = shift;
    my $cache_version = shift;   
    my $use_vep_db = shift;  # 1 for testing/demo, 0 for production
    my $output_is_vep = shift;  # True if vep_output is "vep"

    # assembly and cache_version are optional; if value is empty, not written to config file

    my $output_vep_int = 0;
    if ($output_is_vep) {
        $output_vep_int = 1; 
    }

    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$module.vcf = $vcf
$module.output = $output_fn
$module.vep_cmd = $vep_cmd
$module.cachedir = $cache_dir
$module.reffasta = $reference
$module.usedb = $use_vep_db  
$module.output_vep = $output_vep_int
EOF

    if ($assembly) {
        print OUT "$module.assembly = $assembly\n";
    }
    if ($cache_version) {
        print OUT "$module.cache_version = $cache_version\n";
    }
}

1;
