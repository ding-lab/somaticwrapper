# Annotate VCF file using vcf and write output as VCF or VEP
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
# * preserve_cache_gz

# * vep_output: vcf, vep
# * af_gnomad: path to gnomAD database.  Passed to vep as --af_gnomad.  
# * af_exac: path to ExAC database.  Passed to vep as --af_exac
#
# Annotation is performed using vep.  We can use either a local cache or online
# ('db') lookups (the latter implemented only for vep_output='vcf' or 'vep')..
#
#    VEP Cache logic:
#    * If cache_dir is defined, it indicates location of VEP cache 
#    * if cache_dir is not defined, and cache_gz is defined, extract cache_gz contents into "./vep-cache" and use VEP cache
#    * if neither cache_dir nor cache_gz defined, will perform online VEP DB lookups
#    Online VEP database lookups a) uses online database (so cache isn't installed) b) does not use tmp files
#    It is meant to be used for testing and lightweight applications.  Use the cache for better performance.
#    See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
#
# if extracting cache_gz, then copy it to $cache_dir="./vep-cache" and extract it there
#   (cache_gz is a .tar.gz version of VEP cache; this is typically used for a cwl setup where arbitrary paths are not accessible)
#   These contents will subsequently be deleted, unless preserve_cache_gz is true
#   Cache installation is done in somaticwrapper/image.setup/D_VEP
#   See https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
#
# assembly is the assembly argument passed to vep.  Optional
# vep_output: Defines output format after annotation.  Allowed values: vcf, vep, maf.  Default is vcf
# Additional annotation can be specified with af_exac and af_gnomad arguments, which are passed to vep directly
#   https://useast.ensembl.org/info/docs/tools/vep/script/vep_example.html#gnomad
# 
# Note on multiple annotation, from http://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick
# VEP documentation writes,
#       By default VEP is configured to provide annotation on every genomic feature
#       that each input variant overlaps. This means that if a variant overlaps a
#       gene with multiple alternate splicing variants (transcripts), then a block of
#       annotation for each of these transcripts is reported in the output. In the
#       default VEP output format each of these blocks is written on a single line of
#       output; in VCF output format the blocks are separated by commas in the INFO
#       field.
#       Wherever possible we would discourage users from summarising data in this
#       way. Summarising inevitably involves data loss, and invariably at some point
#       this will lead to the loss of biologically relevant information. For example,
#       if your variant overlaps both a regulatory feature and a transcript and you
#       use one of the flags below, the overlap with the regulatory feature will be
#       lost in your output, when in some cases this may be a clue to the "real"
#       functional effect of your variant. For these reasons we encourage users to
#       use one of the flagging options (--flag_pick, --flag_pick_allele or
#       --flag_pick_allele_gene) and to   post-filter results.
# As a result, we add the --flag_pick flag to vep_opts, then post-filter according to the PICK field
# 
# Output is $results_dir/vep/output.vcf

my $perl = "/usr/bin/perl";  # this is for docker environment

sub vep_annotate {
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $reference = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $assembly = shift;
    my $cache_version = shift; # e.g., 90
    my $cache_dir = shift;  
    my $cache_gz = shift;  
    my $preserve_cache_gz = shift;  # NEW: don't delete cache at end of step if true
    my $vep_output = shift;  # Output format following vep annotation.  May be 'vcf', 'vep', 'maf'
    my $input_vcf = shift;  # Name of input VCF to process
    my $af_exac = shift;  
    my $af_gnomad = shift;  

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

    # if $cache_dir is defined, confirm it exists
    # else if $cache_dir is a .tar.gz file, extract its contents to ./vep-cache
    # else, use VEP DB mode
    if (defined $cache_dir) {
        die "\nError: Cache dir $cache_dir does not exist\n" if (! -d $cache_dir);
        $use_vep_db = 0;
    } elsif ( defined $cache_gz and $cache_gz =~ /\.tar\.gz/ ) {
        $cache_dir = "./vep-cache";
        if (! -d $cache_dir) {
            mkdir $cache_dir or die "$!\n";
        }
        print STDERR "Extracting VEP Cache tarball $cache_gz into $cache_dir\n";
        # This is a preferred way to make system calls - check return value and raise error if necessary
        my $rc = system ("tar -zxf $cache_gz --directory $cache_dir");
        die("Exiting ($rc).\n") if $rc != 0;
        $use_vep_db = 0;
    } else {
        print STDERR "Using online VEP DB\n";
    }

    if ($vep_output =~ /vcf/) {
        $output_fn = "$filter_results/output.vcf";
    } else {
        $output_fn = "$filter_results/output.vep";
    }

    # annotating merged output with VEP annotation
    # followed by AF (allele frequency) and Consequence filter
    write_vep_input(
        $config_fn,          # Config fn
        "merged.vep",                               # Module
        $input_vcf,                # VCF (input)
        $output_fn,           # output
        $vep_cmd, $cache_dir, $reference, $assembly, $cache_version, $use_vep_db, $vep_output =~ /vep/ ? 1 : 0, $af_exac, $af_gnomad);

    my $cmd = <<"EOF1";
    export JAVA_OPTS=\"-Xms256m -Xmx512m\"
    $perl $gvip_dir/vep_annotator.pl $config_fn
EOF1

    my $out = "$job_files_dir/$current_job_file";
    print STDERR "Writing to $out\n";
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
    print STDERR "Executing:\n $bsub_com \n";

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;

    # Clean up by deleting contents of cache_dir - this tends to be big (>10Gb) and unnecessary to keep
    if ( $cache_gz ) {
        if ( $preserve_cache_gz ) {
            print STDERR "Not deleting $cache_dir \n";
        } else {
            print STDERR "Deleting $cache_dir\n";
            my $rc = system("rm -rf $cache_dir\n");
            die("Exiting ($rc).\n") if $rc != 0;
        }
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
    my $af_exac = shift;  
    my $af_gnomad = shift;  

    # assembly and cache_version are optional; if value is empty, not written to config file

    my $output_vep_int = 0;
    if ($output_is_vep) {
        $output_vep_int = 1; 
    }

    print STDERR "Writing to $config_fn\n";
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$module.vcf = $vcf
$module.output = $output_fn
$module.vep_cmd = $vep_cmd
$module.cachedir = $cache_dir
$module.reffasta = $reference
$module.usedb = $use_vep_db  
$module.output_vep = $output_vep_int
$module.vep_opts = --flag_pick
EOF

    # optional parameters
    print OUT "$module.assembly = $assembly\n" if ($assembly);
    print OUT "$module.cache_version = $cache_version\n" if ($cache_version);
    print OUT "$module.af_exac = $af_exac\n" if ($af_exac);
    print OUT "$module.af_gnomad = $af_gnomad\n" if ($af_gnomad);

}

1;
