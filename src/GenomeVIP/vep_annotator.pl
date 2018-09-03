#!/usr/bin/env perl 
#----------------------------------
# @name GenomeVIP VEP annotator script
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# @version 0.2: allow VEP options to be passed
# @version 0.1: original version
#--------------------------------------
use strict;
use warnings;

use Cwd;
use Carp;
use FileHandle;
use IO::File;
use Getopt::Long;
use POSIX qw( WIFEXITED );
use File::Temp qw/ tempfile /;


# if parameter flag output_vep is set, output will be in vep rather than VCF format

# get paras from config file
my (%paras);
map {chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/^\s+//;  $_[1] =~ s/\s+$//; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
map { print STDERR; print STDERR "\t"; print STDERR $paras{$_}; print STDERR "\n" } keys %paras;

# check if options are present
my $opts = "";
my $vcf_flag="--vcf";
if( exists($paras{'vep_opts'}) ) { $opts = $paras{'vep_opts'} };
if( exists($paras{'output_vep'}) && ($paras{'output_vep'}) ) { $vcf_flag = "--symbol" };

# --assembly is optional.  For Cache mode, --cache_version also optional
if (exists($paras{'assembly'})) { $opts = "$opts --assembly $paras{'assembly'}" }

# db mode 1) uses online database (so cache isn't installed) 2) does not use tmp files
# It is meant to be used for testing and lightweight applications.  Use the cache for
# better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
my $cmd="";

if ($paras{'usedb'}) {
    print STDERR ("VEP DB Query mode\n");
    print STDERR ("Reading: $paras{'vcf'}\n");
    print STDERR ("Writing: $paras{'output'}\n");

    # Cannot use --max_af and other AF-related args here
 
    $cmd = "perl $paras{'vep_cmd'} $opts --database --port 3337 --buffer_size 10000  --fork 4 --format vcf $vcf_flag -i $paras{'vcf'} -o $paras{'output'} --force_overwrite  --fasta $paras{'reffasta'}";

    print STDERR $cmd . "\n";
    my $errcode = system($cmd); 
    die ("Error executing: $cmd \n $! \n") if ($errcode);

} else {
    print STDERR ("VEP Cache mode\n");
    if (exists($paras{'cache_version'})) { $opts = "$opts --cache_version $paras{'cache_version'}" }
    my $opts = "$opts --af --max_af --af_1kg --af_esp --af_gnomad";
    $cmd = "perl $paras{'vep_cmd'} $opts --buffer_size 10000 --offline --cache --dir $paras{'cachedir'} --fork 4 --format vcf $vcf_flag -i $paras{'vcf'} -o $paras{'output'} --force_overwrite  --fasta $paras{'reffasta'}";
    print STDERR $cmd . "\n";
    my $errcode = system($cmd); 
    die ("Error executing: $cmd \n $! \n") if ($errcode);
}

print STDERR ("Written to: $paras{'output'}\n");

1;

