#!/usr/bin/perl

## tumor >= 5% and normal <=1% 
### mutect1.7 filtering ###

use strict;
use warnings;
die unless @ARGV == 2;
my ($f_m,$f_filter_out)=@ARGV;

#my $f_m=$run_dir."/mutect/mutect.raw.vcf";
#my $f_filter_out=$run_dir."/mutect/mutect.filtered.vcf";

my $min_vaf_somatic=0.05;
my $max_vaf_germline=0.02;
my $min_coverage=20;

open(IN,"<$f_m");
open(OUT,">$f_filter_out");

while(<IN>) 
	{
		my $l=$_; 
		chomp($l); 
		if($l=~/^#/) { print OUT $l,"\n"; }
		else { 
		my @temp=split("\t",$l); 
		my $tumor=$temp[10]; 
		my $normal=$temp[9]; 
		my @tempt=split(":",$tumor); 
		my @tempn=split(":",$normal); 
		my @readt=split(",",$tempt[1]); 
		my @readn=split(",",$tempn[1]);
		#print $l,"\n";
		#print $readn[0], "\t", $readn[1],"\n";
		#<STDIN>;
		#print $readt[0], "\t", $readt[1],"\n";

		my $tot_n=$readn[0]+$readn[1]; 
		my $tot_t=$readt[0]+$readt[1];										

		if($tot_n==0 || $tot_t==0) { next; }

		my $vaf_t=$readt[1]/$tot_t;
		my $vaf_n=$readn[1]/$tot_n;
		#print $vaf_t,"\t",$vaf_n,"\n";

		if($temp[6] eq "PASS" && $vaf_t>=$min_vaf_somatic && $vaf_n<=$max_vaf_germline && $tot_n>=$min_coverage && $tot_t>=$min_coverage) 
			{
				print OUT $l,"\n"; 
			}	
		}	
	
	}


