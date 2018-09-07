#!/usr/bin/perl

use strict;
use warnings;
die unless @ARGV == 1;

## merge calls from different chromosomes ##

my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
my ($sample_full_path)=@ARGV;
my $f_snv;
my $f_indel; 

my $f_snv_idx; 
my $f_indel_idx; 

my $f_raw; 
my $f_raw_idx; 

my $f_snv_out=$sample_full_path."/mutect1/mutect.filter.snv.vcf";
my $f_ind_out=$sample_full_path."/mutect1/mutect.filter.indel.vcf";

open(OUT1,">$f_snv_out"); 
open(OUT2,">$f_ind_out"); 

foreach my $chr (@chrlist)
    {

	$f_snv=$sample_full_path."/mutect1/mutect.filter.snv.$chr.vcf";
	$f_indel=$sample_full_path."/mutect1/mutect.filter.indel.$chr.vcf";
 
	$f_snv_idx=$sample_full_path."/mutect1/mutect.filter.snv.$chr.vcf.idx";
	$f_indel_idx=$sample_full_path."/mutect1/mutect.filter.indel.$chr.vcf.idx";;

	$f_raw=$sample_full_path."/mutect1/mutect.raw.filtered.$chr.vcf"; 
	$f_raw_idx=$sample_full_path."/mutect1/mutect.raw.filtered.$chr.vcf.idx";
	#$f_raw_idx=$sample_full_path."/mutect/".$sample_name.".gvip.$chr.vcf.idx";	
	#$f_raw=$sample_full_path."/mutect/".$sample_name.".snv.gvip.$chr.vcf";

	if(-e $f_snv) 
	{
		foreach my $l (`cat $f_snv`) 	
		{
		my $ltr=$l; 
		chomp($ltr); 
		if($ltr=~/^#/ &&  !($chr eq "1"))  { next; }
		else { print OUT1 $ltr,"\n"; }
		}
 		`rm $f_snv`;	
		`rm $f_snv_idx`; 
	}

	if(-e $f_indel) 
	{
		foreach my $l (`cat $f_indel`)
    	{
        my $ltr=$l; 
        chomp($ltr);
        if($ltr=~/^#/ &&  !($chr eq "1"))  { next; }
        else { print OUT2 $ltr,"\n"; }
    	}   	
		`rm $f_indel`; 
		`rm $f_indel_idx`; 
 	}

	`rm $f_raw`; 
	#`rm $f_raw`; 
	#`rm $f_raw_idx`;
	`rm $f_raw_idx`; 

    #print GATK "rawvcf=".$sample_full_path."/mutect/".$sample_name.".raw.$chr.vcf\n";
    #print GATK "gvipvcf=".$sample_full_path."/mutect/".$sample_name.".gvip.$chr.vcf\n";
    #print GATK "snvvcf=".$sample_full_path."/mutect/".$sample_name.".snv.gvip.$chr.vcf\n";
    #print GATK "indelvcf=".$sample_full_path."/mutect/".$sample_name.".indel.gvip.$chr.vcf\n";
    #print GATK "     ".$run_script_path."genomevip_label.pl GATK \${rawvcf} \${gvipvcf}"."\n";
	}
