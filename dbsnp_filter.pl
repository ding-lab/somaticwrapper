#!/usr/bin/env perl
#----------------------------------
# @name GenomeVIP dbSNP annotation and filtering script
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
#
# @version 0.3 (rjm): add mode switch
# @version 0.2 (rjm): workaround for stat on AWS
# @version 0.1 (rjm): based on approach from Venkata Yellapantula
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
#use File::stat;


sub checksize {
    my ($dest,$src) = @_;
    my $fs = `wc -l < $dest`;
    if ( $fs == 0){ system("grep ^# $src > $dest"); }
    return 1;
}

# get paras from config file
my (%paras);
map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; print $v."\n";  $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
# map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;

# Use uncompressed db to avoid being bitten by java compression bug
my $anno=$paras{'rawvcf'}."dbsnp_anno.vcf";
if ($paras{'rawvcf'} =~ /\.vcf$/) {
    ($anno = $paras{'rawvcf'}) =~ s/\.vcf$/\.dbsnp_anno\.vcf/;
}


my $cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} annotate -id $paras{'db'} $paras{'rawvcf'} > $anno";
print "$cmd\n";
system($cmd);
checksize($anno, $paras{'rawvcf'});
if( exists $paras{'mode'}  &&  $paras{'mode'} eq "filter" )  {
$cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} filter -n \" (exists ID) & (ID =~ 'rs' ) \" -f $anno > $paras{'passfile'}";
system($cmd);
checksize($paras{'passfile'}, $anno);
$cmd = "java $ENV{'JAVA_OPTS'} -jar $paras{'annotator'} filter    \" (exists ID) & (ID =~ 'rs' ) \" -f $anno > $paras{'dbsnpfile'}";
system($cmd);
checksize($paras{'dbsnpfile'}, $anno);
$cmd = "rm -f $anno";
system($cmd);
}
1;

