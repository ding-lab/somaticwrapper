######### Song Cao###########
#  qc_vcf.pl #

#!/usr/bin/perl
use strict;
use warnings;

die unless @ARGV == 2;

my ($vcf_in,$vcf_out)=@ARGV; 

open(IN,"<$vcf_in"); 
open(OUT,">$vcf_out"); 
my @t; 
my $ref; 
my $var; 
while(<IN>)
{
	my $ltr=$_; 
	chomp($ltr); 
 	if($ltr=~/^#/) 
	{
	print OUT $ltr,"\n"; 	
	}	
	else 
	{
	@t=split("\t",$ltr); 
	# filter case ref eq var ##
	$ref=$t[3]; 
	$var=$t[4];
	if($ref eq $var) 
	{
	next; 
	}
	else { print OUT $ltr,"\n"; }
	}
}

close IN; 
close OUT; 
#my @temp = split("/", $dir);
#my $run_name = pop @temp;
#my $f_out=$run_dir."/Analysis_Summary_Somatic_$run_name.tsv";

