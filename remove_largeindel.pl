#!/usr/bin/perl

use strict;
use warnings;
die unless @ARGV == 2;
my ($file_in,$file_out)=@ARGV;

open(IN,"<$file_in");
open(OUT,">$file_out"); 
my @t;
my $lref; 
my $lvar; 
while(<IN>)
{
	my $l=$_; 
	chomp($l);
	if($l=~/^#/) { print OUT $l, "\n"; next; }
	else {
	@t=split("\t",$l); 
	#$l=~s/SVTYPE=//g;
	$lref=length($t[3]); 
	$lvar=length($t[4]);
	if(abs($lref-$lvar)<=100 && $lref<=101 && $lvar<=101)
	{  
	print OUT $l,"\n"; 
	}
	}
}

close IN; 
close OUT;
