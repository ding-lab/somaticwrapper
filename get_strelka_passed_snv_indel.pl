#!/usr/bin/perl

### 10/10/2018 ##
### get passed snv and indel variants 
### add the filtering for indel length ##

use strict;
use warnings;

die unless @ARGV == 4;
my ($f_snv_in,$f_indel_in,$f_snv_out,$f_indel_out)=@ARGV; 

open(OUT1,">$f_snv_out"); 
open(OUT2,">$f_indel_out");
my $head;
#my $firsttime=1;

foreach my $l (`zcat $f_snv_in`) 
	{
		my $ltr=$l;
		chomp($ltr);  
		if($ltr=~/^#/) { 
#		if($ltr=~/^#CHROM/) 
#		{
#			$head=$ltr; 
#		}
		print OUT1 $ltr,"\n"; 
		next;  }
		else {
		 my @temp=split("\t",$ltr); 
		 #my $ref=$temp[3];
		 #my $var=$temp[4];
		 if($temp[6] eq "PASS")
		 {
			print OUT1 $ltr, "\n"; 
		 }
	}
	}
foreach my $l (`zcat $f_indel_in`) 
    {
        my $ltr=$l;
        chomp($ltr);  
        if($ltr=~/^#/) { print OUT2 $ltr,"\n"; next;  }
        else {
         my @temp=split("\t",$ltr); 
         #my $ref=$temp[3];
         #my $var=$temp[4];
         if($temp[6] eq "PASS")
         {
           #if($firsttime==1) { print OUT2 $head,"\n"; $firsttime=0; }  
		   print OUT2 $ltr, "\n";
         }
   		 }	
	}

close OUT1; 
close OUT2; 

