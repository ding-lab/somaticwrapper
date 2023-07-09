#!/usr/bin/env perl 
#----------------------------------
# @name GenomeVIP VEP annotator script
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# @version 0.2: allow VEP options to be passed
# @version 0.1: original version
#--------------------------------------
## song changed the parameter for vep 
### 12/13/17 ###
use strict;
use warnings;

use Cwd;
use Carp;
use FileHandle;
use IO::File;
use Getopt::Long;
use POSIX qw( WIFEXITED );
use File::Temp qw/ tempfile /;

# get paras from config file
my (%paras);
map {chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/^\s+//;  $_[1] =~ s/\s+$//; my $v = $_[1]; print $v."\n";  $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
 map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;

# check if options are present
my $opts="";
if( exists($paras{'vep_opts'}) ) { $opts = $paras{'vep_opts'} };

my $cmd="";

# split off original header
my (undef, $tmp_orig_calls)  = tempfile();
$cmd="/bin/grep -v ^# $paras{'vcf'} > $tmp_orig_calls";
   system($cmd);

# run vep
my (undef, $tmp_vep_out) = tempfile();

### wxs less mutation ##
#$cmd = "perl $paras{'vep_cmd'}  --offline --cache --dir $paras{'cachedir'}  --polyphen b --gmaf --maf_1kg --maf_esp --assembly $paras{'assembly'} --fork 4 --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --format vcf --vcf -i $tmp_orig_calls -o $tmp_vep_out --force_overwrite  --fasta $paras{'reffasta'}";

#$cmd = "perl $paras{'vep_cmd'}  --offline --cache --dir $paras{'cachedir'}  --polyphen b --gmaf --maf_1kg --maf_esp --assembly $paras{'assembly'} --fork 4 --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --format vcf --input_file $paras{'vcf'} --output_file $paras{'output'} --force_overwrite --fasta $paras{'reffasta'}";

$cmd = "/usr/bin/perl $paras{'vep_cmd'} $opts --buffer_size 300 --offline --cache --dir $paras{'cachedir'} --assembly $paras{'assembly'} --fork 8 --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --format vcf --vcf -i $tmp_orig_calls -o $tmp_vep_out --force_overwrite  --fasta $paras{'reffasta'}";

## for wgs simple ##

#$cmd = "perl $paras{'vep_cmd'} $opts --buffer_size 10000 --offline --cache --dir $paras{'cachedir'} --assembly $paras{'assembly'} --fork 4 --format vcf --vcf -i $tmp_orig_calls -o $tmp_vep_out --force_overwrite  --fasta $paras{'reffasta'}";

#$cmd = "perl $paras{'vep_cmd'} $opts --buffer_size 1000 --offline --cache --dir $paras{'cachedir'} --assembly $paras{'assembly'} --fork 4 --format vcf --vcf -i $tmp_orig_calls -o $tmp_vep_out --force_overwrite --everything --fasta $paras{'reffasta'}";

#variant_effect_predictor.pl --species $species --assembly $ncbi_build --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir $vep_data --fasta $ref_fasta --format vcf --input_file $input_vcf --output_file $output_vcf
 
system($cmd);

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

1;

