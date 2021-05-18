#!/usr/bin/perl

## tumor vaf >=5% 
### mutect2 filtering ###

use strict;
use warnings;
die unless @ARGV == 4;

my ($f_filter_in,$f_filter_out,$mincov_t,$minvaf)=@ARGV;


my $min_vaf_somatic=$minvaf;
#my $min_coverage=$mincov;

open(IN,"<$f_filter_in");
open(OUT,">$f_filter_out");

while(<IN>) 
	{
		my $l=$_; 
		chomp($l); 
		#print $l,"\n";
		if($l=~/^#/) { #print OUT $l,"\n";
   		if($l=~/^#CHROM/) {  
		my @temphead=split("\t",$l); 
		print OUT $temphead[0]; 
		for(my $i=1;$i<=8;$i++) 
		{
		print OUT "\t",$temphead[$i]; 
		}
		print OUT "\t","TUMOR","\n"
		}
		else { print OUT $l,"\n"; } 
 		}
		else { 
		my @temp=split("\t",$l); 

		my $tumor=$temp[9];
		#print $tumor,"\n"; <STDIN>;
		my @tempt=split(":",$tumor); 
		#my @tempn=split(":",$normal); 
		my @readt=split(",",$tempt[1]); 

#### readt[0]: read count for ref allele in tumor ##
### readt[1]: read count for alt allele in tumor ##
 		
		#my $tot_n=$readn[0]+$readn[1]; 
		my $tot_t=$readt[0]+$readt[1];										

		if($tot_t==0) { next; }
		
		my $vaf_t=$readt[1]/$tot_t;

#		my $vaf_t=$tempt[2];
		#print $tot_t,"\n";
		#print $vaf_t,"\n";
	## apply filtering ##
		if($temp[6] eq "PASS" &&  $vaf_t>=$minvaf && $tot_t>=$mincov_t) 
			{
				print OUT $l,"\n"; 
			}	
		}	
	
}

close IN;
close OUT; 
