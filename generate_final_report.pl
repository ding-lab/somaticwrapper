#!/usr/bin/perl

### for example, 
## tumor >= 5% and normal <=1% 
use strict;
use warnings;
die unless @ARGV == 2;
my ($run_dir, $s_exonic)=@ARGV;

my $working_name= (split(/\//,$run_dir))[-1];

my $f_sum=$run_dir."/".$working_name.".withmutect.maf\n";
my $f_status=$run_dir."/".$working_name.".withmutect.status\n";

open(OUT1,">$f_sum");
open(OUT2,">$f_status"); 

my $head_w=0;

foreach my $d (`ls $run_dir`)
{
  	my $dtr=$d; 
	chomp($dtr);
 
#	my $f_maf=$run_dir."/".$dtr."/".$dtr.".checked.maf"; 
my $f_maf=$run_dir."/".$dtr."/".$dtr.".withmutect.filtered.maf"; 
	if(-e $f_maf) 
	{
		my $count=0;
		foreach my $l (`cat $f_maf`) 
		{ 
			my $ltr=$l; 
			chomp($ltr); 
			if($ltr=~/^#version/)  { next; } 
			else { 
			if($ltr=~/^Hugo/) { if($head_w==0) { print OUT1 $ltr,"\n"; } $head_w=1; } 
			else { 
				my @temp=split("\t",$ltr); 
				my $annot=$temp[8];
				my $af=$temp[99];
				if($annot=~/Frame_Shift_Del/ || $annot=~/Frame_Shift_Ins/ || $annot=~/Missense_Mutation/ || $annot=~/Nonsense_Mutation/ ||  $annot=~/Nonstop_Mutation/ || $annot=~/Silent/ || $annot=~/Splice_Site/ || $annot=~/In_Frame_Ins/ || $annot=~/In_Frame_Del/  || $s_exonic==0) {
					#if($af eq "" || (($af ne "") && $af<0.005))
					#{
					print OUT1 $ltr,"\n"; 
					$count++;
					#}
					} 
				} 
        		}
      	}
		print OUT2 $count,"\t",$f_maf,"\n";
	} 
 
} ##


close OUT1;
close OUT2; 

