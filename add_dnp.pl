#!/usr/bin/perl


## Song Cao ##

use strict;
use warnings;

die unless @ARGV == 3;

##adding correct annotation for dnp in the maf file ##
##fix dnps which are in different codon and same allele ##

## use mnp for general purpose ##
 
my %mnp_infor=(); 
my %mnp_end=(); 

my ($maf_in,$maf_mnp,$maf_out)=@ARGV;

open(IN,"<$maf_in");
open(MNP,"<$maf_mnp");
open(OUT,">$maf_out"); 

### merged ###

my $status_mnp=0; 
my $ref=""; 
my $var=""; 
my $prot1="";
my $prot2="";
my $effect="";
my $c=""; 

foreach my $l (`cat $maf_mnp`) 
{ 
	my $ltr=$l; 
	chomp($ltr);
	my @l=split("\t",$ltr);
 
	if($ltr=~/^Hugo/) { next; } 
	else 
	{ 

	if($ltr=~/^#/)
	{
		if($ltr=~/Merged\=Yes/) 
		{ $status_mnp=1; } 
		elsif($ltr=~/2\s+1\s+Merged=No\(different codons\)/) 
		{ 
			$status_mnp=2;
			$ref="";
			$var="";
			$prot1="";
			$prot2="";
			$effect="";
			$c="";
			#print $ltr,"\n";
			#<STDIN>;
		 } 
		elsif($ltr=~/2\s+1\s+Merged=No\(not cds\)/)
		{ 
			$status_mnp=2; 
			$ref="";
			$var=""; 
			$prot1="";
			$prot2="";
			$effect="";
			$c="";
			#print $ltr,"\n";
            #<STDIN>;
		}else { 
		$status_mnp=0; 
		}
	}
	else 
	{ 
		if($status_mnp==1)  
		{
			#my @l=split("\t",$ltr); 
			my $sn=$l[15];
			$mnp_infor{$sn}{$l[4]}{$l[5]}=$ltr; 											
			$mnp_end{$sn}{$l[4]}{$l[5]}=$l[6];		
		}

		if($status_mnp==2) 
		{
			$ref.=$l[10]; 
			$var.=$l[12];
			
			if($prot2 eq "")
			{
			$c=$l[34];
			$prot1=$l[35];
			$prot2=$l[36];
			}
			else 
			{
			$c.=";".$l[34];
			$prot1.=";".$l[35];
			$prot2.=";".$l[36];
			}

			#print $l[9],"\n";

            if($l[8]=~/Splice_Site/) { $effect=$l[8]; }
			## silent ##
            if($l[8]=~/Silent/ && !($effect=~/Nonsense_Mutation/ || $effect=~/Missense_Mutation/ || $effect=~/Splice_Site/ || $effect=~/Nonstop_Mutation/))                 
			{
            	$effect=$l[8];
            }
   			
            if($l[8]=~/Missense_Mutation/ && !($effect=~/Nonsense_Mutation/ || $effect=~/Splice_Site/ || $effect=~/Nonstop_Mutation/)) {
                $effect=$l[8];
            }

		    if($l[8]=~/Nonstop_Mutation/ && !($effect=~/Nonsense_Mutation/ || $effect=~/Splice_Site/)) {
                $effect=$l[8];
            }	

		    if($l[8]=~/Nonsense_Mutation/ && !($effect=~/Splice_Site/)) {
                $effect=$l[8];
            } 

			if(length($ref)==2 && length($var)==2)
			{
			my $start=$l[5]-1; 
			my $end=$l[5];			

			#print $ref,"\t",$var,"\t",$effect,"\t",$start,"\t",$c,"\t",$prot1,"\t",$prot2,"\n";
			#<STDIN>;
			$mnp_infor{$l[15]}{$l[4]}{$start}=join("\t", @l[0..4])."\t$start\t$end\t$l[7]\t$effect\tMNP\t$ref\t$var\t$var\t".join("\t",@l[13..33])."\t".$c."\t".$prot1."\t".$prot2."\t".join("\t", @l[37..$#l]);
			$mnp_end{$l[15]}{$l[4]}{$start}=$end;
			} 
				  				
		}
	
	}
 
	} 
}

my $removed=0; 
my $chr;
my $begin;  
my $end; 

foreach my $l (`cat $maf_in`) 
	{

		my $ltr=$l; 
		chomp($ltr); 
		if($ltr=~/^Hugo/) { print OUT $ltr,"\n"; }
		else 
		{ 
		my @l=split("\t",$ltr); 
		if($ltr=~/^Hugo/) { print OUT $ltr,"\n"; }
		else { 
		my @l=split("\t",$ltr);
		my $sn=$l[15]; 	

		## mnp ##

		if(defined $mnp_infor{$sn}{$l[4]}{$l[5]}) 
		{
			print OUT $mnp_infor{$sn}{$l[4]}{$l[5]},"\n";
			#print $sn,"\t",$l[4],"\t",$l[5],"\n"; 
			#print $mnp_infor{$sn}{$l[4]}{$l[5]},"\n";
			#<STDIN>;
		 	$removed=1; 
			$chr=$l[4]; 
			$begin=$l[5]; 
			$end=$mnp_end{$sn}{$l[4]}{$l[5]}									
		}
		else 
		{ 
		if($removed==1 && $l[4] eq $chr && ($l[5]>=$begin && $l[5]<=$end)) 
		{
		#print $l[4],"\t",$l[5],"\n";
		#<STDIN>;
		next; 				
		}
		else { print OUT $ltr,"\n"; $removed=0; }
		}
																		
		}		
		}

	}


close IN; 
close MNP;
close OUT; 
