#!/usr/bin/perl

## tumor >= 5% and normal <=1% 
### add the filtering for indel length ##
### mutect filtering ###

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
		my $tumor=$temp[9]; 
		my $normal=$temp[10]; 
		my @tempt=split(":",$tumor); 
		my @tempn=split(":",$normal); 
		my @readt=split(",",$tempt[1]); 
		my @readn=split(",",$tempn[1]);
		my $tot_n=$readn[0]+$readn[1]; 
		my $tot_t=$readt[0]+$readt[1];										
		my $vaf_t=$tempt[2];
		my $vaf_n=$tempn[2];
		if($temp[6] eq "PASS" && $vaf_t>=$min_vaf_somatic && $vaf_n<=$max_vaf_germline && $tot_n>=$min_coverage && $tot_t>=$min_coverage) 
			{
				print OUT $l,"\n"; 
			}	
		}	
	
	}


