#!/usr/bin/perl

use warnings;
die unless @ARGV == 2;
my ($f_in,$f_out)=@ARGV;

#my $f_m=$run_dir."/mutect/mutect.raw.vcf";
#my $f_filter_out=$run_dir."/mutect/mutect.filtered.vcf";

open(IN,"<$f_in");
open(OUT,">$f_out");

while(<IN>) 
	{
		my $l=$_; 
		chomp($l);
		$l=~s/\tNaN\t1\.0/\t0\t0/;  
		print OUT $l,"\n"; 	
	}
close IN; 
close OUT;
