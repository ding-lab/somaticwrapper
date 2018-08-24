# Annotate VCF file and write output as MAF, using vcf2maf.pl
#
# Required Arguments (from command line):
# * input_vcf 
# * reference_fasta 
# 
# Optional Arguments
# * assembly
# * cache_version
# * cache_dir
# * cache_gz
# * exac_vcf: path to ExAC database.  Passed to vcf2maf.pl as --filter-vcf. 0 to disable

# Annotation is performed using [vcf2maf](https://github.com/mskcc/vcf2maf), which uses vep.
#    --cache_gz: is a file ending in .tar.gz containing VEP cache tarball
#    --cache_dir s: location of VEP cache directory
#        VEP Cache logic:
#        * If cache_dir is defined, it indicates location of VEP cache 
#        * if cache_dir is not defined, and cache_gz is defined, extract cache_gz contents into "./vep-cache" and use VEP cache
#        * if neither cache_dir nor cache_gz defined, error.  vcf_2_maf does not support online cache lookups

# if $cache_gz file exists, then copy it to $cache_dir="./vep-cache" and extract it there
#   (assuming it is a .tar.gz version of VEP cache; this is typically used for a cwl setup where arbitrary paths are not accessible)
#   These contents will subsequently be deleted
#   Cache installation is done in somaticwrapper/image.setup/D_VEP
#   See https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
#
# assembly is passed as --ncbi-build to vcf2maf.pl
# Additional annotation can be specified with --exac_vcf argument, which is passed as `--filter-vcf` argument to vcf2maf.pl
#   Pass value of 0 to disable.
# see vcf_2_maf.md for additional details
# 
# Output is $results_dir/maf/output.maf


my $perl = "/usr/bin/perl";  # this is for docker environment

sub vcf_2_maf {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $reference = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    # assembly and cache_version may be blank; if so, not passed on command line to vep
    my $assembly = shift;
    my $cache_version = shift; # e.g., 90
    my $cache_dir = shift;  
    my $cache_gz = shift;  
    my $input_vcf = shift;  # Name of input VCF to process
    my $exac_vcf = shift;   # passed as --filter-vcf

    $current_job_file = "j9_vcf_2_maf.sh";

    my $bsub = "bash";
    my $filter_results = "$results_dir/maf";
    system("mkdir -p $filter_results");

    $output_fn = "$filter_results/output.maf";

    # if $cache_dir is defined, confirm it exists
    # if $cache_dir is not defined, and $cache_gz is a .tar.gz file, extract its contents to ./vep-cache
    if (defined $cache_dir) {
        die "\nError: Cache dir $cache_dir does not exist\n" if (! -d $cache_dir);
    } elsif ( defined $cache_gz and $cache_gz =~ /\.tar\.gz/ ) {
        $cache_dir = "./vep-cache";
        if (! -d $cache_dir) {
            mkdir $cache_dir or die "$!\n";
        }
        print STDERR "Extracting VEP Cache tarball $cache_gz into $cache_dir\n";
        # This is a preferred way to make system calls - check return value and raise error if necessary
        my $rc = system ("tar -zxf $cache_gz --directory $cache_dir");
        die("Exiting ($rc).\n") if $rc != 0;
    } else {
        die "--cache_dir or --cache_gz must be defined\n";
    }


#      --ncbi-build     NCBI reference assembly of variants MAF (e.g. GRCm38 for mouse) [GRCh37]
#      --cache-version  Version of offline cache to use with VEP (e.g. 75, 84, 91) [Default: Installed version]
#      --vep-data       VEP's base cache/plugin directory [~/.vep]

    my $opts = "--vep-data $cache_dir";
    if ($assembly) { $opts = "$opts --ncbi-build $assembly"; }
    if ($cache_version)  { $opts = "$opts --cache-version $cache_version"; }

    my $vep_path = dirname($vep_cmd);
    my $cmd = "$perl /usr/local/mskcc-vcf2maf/vcf2maf.pl $opts --input-vcf $input_vcf --output-maf $output_fn --ref-fasta $reference --filter-vcf $exac_vcf --vep-path $vep_path --tmp-dir $filter_results";

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


1;
