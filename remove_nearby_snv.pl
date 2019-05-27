#!/usr/bin/perl


## Song Cao ##

use strict;
use warnings;

die unless @ARGV == 2;

## chr: 4, pos 5, type 9, sn 15 starting from 0 ##

###  762 DEL
### 167 INS
## 7938 SNP ###

my %ind=(); 
my %snv=(); 
my $window=20; 

my ($maf_in,$maf_out)=@ARGV;

open(OUT,">$maf_out"); 
my $maf_out_rem=$maf_out.".removed"; 
open(OUT2,">$maf_out_rem"); 

foreach my $l (`cat $maf_in`) 
{ 
	my $ltr=$l; 
	chomp($ltr); 
	if($ltr=~/^Hugo/) { next; } 
	else { 
	my @temp=split("\t",$ltr);
	my $id=$temp[15]."_".$temp[4]; 

	if($temp[9] eq "INS" || $temp[9] eq "DEL") 
	{
		$ind{$id}{$temp[5]}=1;		
	}

	if($temp[9] eq "SNP")
        {
        	$snv{$id}{$temp[5]}=1;   
        }
	
	} 
}

foreach my $l (`cat $maf_in`) 
	{
		my $ltr=$l; 
		chomp($ltr); 
		if($ltr=~/^Hugo/) { print OUT $ltr,"\n"; }
		else { 
		my @temp=split("\t",$ltr); 
		my $id=$temp[15]."_".$temp[4];	
		if($temp[9] eq "SNP") 
		{
			my $remove_snv=0; 
			if(defined $ind{$id}) 
			{
			for(my $i=-$window;$i<=$window;$i++)
			{
			my $pos=$temp[5]+$i; 
			if(defined $ind{$id}{$pos}) 
			{
			  $remove_snv=1; 		
			}			
			}
			if($remove_snv==0) 
			{
			print OUT $ltr,"\n"; 	
			}
			if($remove_snv==1) 
			{
			print OUT2 $ltr,"\n"; 
			}	
			}
			else { print OUT $ltr,"\n"; }	
		}
		else { print OUT $ltr,"\n"; }			
		} 
	}

close OUT; 
close OUT2; 
