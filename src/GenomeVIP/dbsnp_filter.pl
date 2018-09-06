#!/usr/bin/env perl
#----------------------------------
# @name GenomeVIP dbSNP annotation and filtering script
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# @author Matthew A. Wyczalkowski  <m.wyczalkowski@genome.wustl.edu>
#
# @version 0.4 (maw): documentation updates
# @version 0.3 (rjm): add mode switch
# @version 0.2 (rjm): workaround for stat on AWS
# @version 0.1 (rjm): based on approach from Venkata Yellapantula
#--------------------------------------
use strict;
use warnings;

#use Cwd;
#use Carp;
#use FileHandle;
#use IO::File;
#use Getopt::Long;
#use POSIX qw( WIFEXITED );
#use File::Temp qw/ tempfile /;
#use File::stat;
use File::Basename;

# paras{} arguments used:
# * rawvcf - the input VCF file
# * annotator - snpSift JAR file
# * db - dbsnp database file (vcf)
# * mode - should be 'filter', else annotate without filtering
# * passfile - output file with entries NOT in dbsnp_db
# * presentfile -  output file with entries IN dbsnp_db

# This function seems to add header lines from src to dest if dest is empty.
sub checksize {
    my ($dest,$src) = @_;
    my $fs = `wc -l < $dest`;
    if ( $fs == 0){ system("grep ^# $src > $dest"); }
    return 1;
}

# get paras from config file
my (%paras);
map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
map { print STDERR; print STDERR "\t"; print STDERR $paras{$_}; print STDERR "\n" } keys %paras;

# 1. use SnpSift to add info to "ID" column based on whether variant exists in 'db'
#    implicitly assume all variant names begin with 'rs'
#    This is dbsnp_anno file
# if mode == filter:
#   2. Create 'passfile' with all entries which do NOT have ID field 
#   3. Create 'presentfile' with all entries which DO have ID field
# else
#   a. Create 'passfile' as link to dbsnp_anno

# Create annotation file; this is a temporary file unless not filtering.
# Annotation file will be in same directory as dbsnp_pass

my($filename_in, $dirs_in, $suffix_in) = fileparse($paras{'rawvcf'});
my($filename_pass, $dirs_pass, $suffix_pass) = fileparse($paras{'passfile'});
my $anno = $dirs_pass . "/" . $filename_in . ".dbsnp_anno.vcf";

my $cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} annotate -id $paras{'db'} $paras{'rawvcf'} > $anno";
print STDERR "$cmd\n";
my $return_code = system ( $cmd );
die("Exiting ($return_code).\n") if $return_code != 0;
checksize($anno, $paras{'rawvcf'});

if( exists $paras{'mode'}  &&  $paras{'mode'} eq "filter" )  {
    $cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} filter -n \" (exists ID) & (ID =~ 'rs' ) \" -f $anno > $paras{'passfile'}";
    print STDERR "$cmd\n";
    my $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;
    checksize($paras{'passfile'}, $anno);

    $cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} filter    \" (exists ID) & (ID =~ 'rs' ) \" -f $anno > $paras{'presentfile'}";
    print STDERR "$cmd\n";
    $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;
    checksize($paras{'presentfile'}, $anno);

    print STDERR "Deleting temporary file $anno\n";
    $cmd = "rm -f $anno";
    $return_code = system ( $cmd );
    die("Exiting ($return_code).\n") if $return_code != 0;
} else {
    my $anno_abs = `readlink -f $anno`;
    chomp $anno_abs;
    print STDERR "Making link to $anno_abs from $paras{'passfile'}\n";
    $return_code = system ("ln -fs $anno_abs $paras{'passfile'} ");
    die("Exiting ($return_code).\n") if $return_code != 0;
}
1;

