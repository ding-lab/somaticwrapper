#!/usr/bin/env perl
# ---------------------
# @name GenomeVIP utility to extract 'other' calls from varscan
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# @version
#
# @syntax  extract_somatic_other.pl  < source.vcf  >  destination.vcf
# ---------------------
use strict;
use warnings;


while(<>) {
    my @a = split /\t/;
    print $_  if( /^#/ || ($a[7] =~ /^SS=[05]/ ||  $a[7] =~ /;SS=[05]$/ ||  $a[7] =~ /;SS=[05];/) );
}
