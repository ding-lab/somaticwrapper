#!/usr/bin/perl

### for example, 
use strict;
use warnings;
die unless @ARGV == 2;
my ($f_in, $f_out)=@ARGV;

open(OUT,">$f_out"); 
my %count_v; 
my $f_status=$f_out.".status\n";

open(OUT1,">$f_status"); 

foreach my $l (`cat $f_in`) 
    { 
	my $ltr=$l; 
	chomp($ltr); 
	if($ltr=~/^#version/)  { next; } 
	else { 
		if($ltr=~/^Hugo/) { print OUT $ltr,"\n"; next; } 
			else { 
				my @temp=split("\t",$ltr); 
				my $annot=$temp[8];
				my $sn=$temp[15]; 
				$sn=~s/_T//g; 
				if($annot=~/Frame_Shift_Del/ || $annot=~/Frame_Shift_Ins/ || $annot=~/Missense_Mutation/ || $annot=~/Nonsense_Mutation/ ||  $annot=~/Nonstop_Mutation/ || $annot=~/Silent/ || $annot=~/Splice_Site/ || $annot=~/In_Frame_Ins/ || $annot=~/In_Frame_Del/) {
					print OUT $ltr,"\n"; 
					$count_v{$sn}++; 
					} 
				} 
        		}
      	}

close OUT;

foreach my $s (sort keys %count_v)
{
	print OUT1 $s,"\t",$count_v{$s},"\n"; 
}

close OUT1; 
