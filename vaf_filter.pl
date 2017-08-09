#!/usr/bin/perl

## tumor >= 5% and normal <=1% 
use strict;
use warnings;
die unless @ARGV == 1;
my ($run_dir)=@ARGV; 

my $f_m=$run_dir."/merged.vcf"; 
my $f_filter_out=$run_dir."/merged.filtered.vcf";
my $f_vaf_out=$run_dir."/merged.vaf";

open(OUT1,">$f_filter_out");
open(OUT2,">$f_vaf_out"); 

foreach my $l (`cat $f_m`) 
	{
		my $ltr=$l; 
		if($ltr=~/^#/) { next;  }
		else {
	 	my $info=$temp[7];
		#my $r_info=$temp[8]; 
		 
		if($info=~/strelka/) 
		{
		my $vaf_n=$temp[11];
		my $vaf_t=$temp[12];
		my $ref=$temp[3]; 
		my $var=$temp[4];
		my $r_tot=0; 	
		my @temp2=split(":",$vaf_n);
		my %rc=();			
		
	 	 $rc{'A'}=split(",",$temp2[0])[0]; 
		 $rc{'C'}=split(",",$temp2[1])[0];
		 $rc{'G'}=split(",",$temp2[4])[0];
			$rc{'T'}=split(",",$temp2[7])[0];

		foreach my $nt (keys %rc) 
		{
			$r_tot+=$rc{$nt}; 
		}

		@temp2=split(":",$vaf_t);

        my %rc2=();
		my $r_tot2=0;

         $rc2{'A'}=split(",",$temp2[0])[0];
         $rc2{'C'}=split(",",$temp2[1])[0];
         $rc2{'G'}=split(",",$temp2[4])[0];
         $rc2{'T'}=split(",",$temp2[7])[0];

        foreach my $nt (keys %rc)
        {
            $r_tot2+=$rc{$nt};
        }
		
		print $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t",$rc{$ref},"\t",$rc{$ref}/$r_tot,"\t",$rc{$var},"\t",$rc{$var}/$r_tot,"\t",$rc2{$ref},"\t",$rc2{$ref}/$r_tot2,"\t",$rc2{$var},"\t",$rc2{$var}/$r_tot2,"\n"; 

		}
	
		if($info=~/varscan/
		{
			
		}
	}
