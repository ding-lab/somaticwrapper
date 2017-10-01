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
map {chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/^\s+//;  $_[1] =~ s/\s+$//; my $v = $_[1]; print $v."\n";  $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
 map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;

# check if options are present
my $opts = "";
my $vcf_flag="--vcf";
if( exists($paras{'vep_opts'}) ) { $opts = $paras{'vep_opts'} };
if( exists($paras{'output_vep'}) && ($paras{'output_vep'}) ) { $vcf_flag = "--symbol" };

# db mode 1) uses online database (so cache isn't installed) 2) does not use tmp files
# It is meant to be used for testing and lightweight applications.  Use the cache for
# better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
my $cmd="";

if ($paras{'usedb'}) {
    print("VEP DB Query mode\n");
    print("Reading: $paras{'vcf'}\n");
    print("Writing: $paras{'output'}\n");
 
    $cmd = "perl $paras{'vep_cmd'} $opts --database --port 3337 --buffer_size 10000  --assembly $paras{'assembly'} --fork 4 --format vcf $vcf_flag -i $paras{'vcf'} -o $paras{'output'} --force_overwrite  --fasta $paras{'reffasta'}";

    print($cmd . "\n");
    system($cmd);

} else {
    print("VEP Cache mode\n");

    # split off original header
    my (undef, $tmp_orig_calls)  = tempfile();
    $cmd="/bin/grep -v ^# $paras{'vcf'} > $tmp_orig_calls";
       system($cmd);

    # run vep if input VCF not empty
    my (undef, $tmp_vep_out) = tempfile();
    if (-s $tmp_orig_calls) {
        $cmd = "perl $paras{'vep_cmd'} $opts --buffer_size 10000 --offline --cache --dir $paras{'cachedir'} --assembly $paras{'assembly'} --fork 4 --format vcf $vcf_flag -i $tmp_orig_calls -o $tmp_vep_out --force_overwrite  --fasta $paras{'reffasta'}";
        print($cmd . "\n");
        system($cmd);
    } else {
        print("VCF is empty\n");
        system("touch $tmp_vep_out");
    }

    # re-merge headers and move
    my (undef, $tmp_merge) = tempfile();
    $cmd = "grep ^##fileformat $tmp_vep_out > $tmp_merge";
       system($cmd);
    $cmd = "grep ^# $paras{'vcf'} | grep -v ^##fileformat | grep -v ^#CHROM >> $tmp_merge";
       system($cmd);
    $cmd = "grep -v ^##fileformat $tmp_vep_out >> $tmp_merge";
       system($cmd);
    $cmd = "cat $tmp_merge > $paras{'output'}";
       system($cmd);

    #Save other output
    my @suffix=("_summary.html", "_warnings.txt");
    foreach (@suffix) {
        my $file = $tmp_vep_out.$_;   if ( -e $file ) { $cmd = "cat $file > $paras{'output'}".$_; system($cmd); }
    }
    # clean up
    $cmd = "rm -f $tmp_orig_calls $tmp_vep_out"."*"." ".$tmp_merge;
    system($cmd);
}

1;

