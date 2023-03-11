#!/usr/bin/perl
### check splice_site if it is correctly annotated ### 
use strict;
use warnings;
die unless @ARGV == 1;
my ($run_dir)=@ARGV;

my $sn=(split(/\//,$run_dir))[-1];

my %splice_s=(); 

my $f_bed="/gscmnt/gc2524/dinglab/bed_maker/E75_bed_v3.tsv"; 
my $f_in=$run_dir."/".$sn.".maf";
my $f_out=$run_dir."/".$sn.".checked.maf";

open(IN,"<$f_bed");
open(OUT,">$f_out"); 

while(<IN>)
{
	my $l=$_; 
	chomp($l);
	my @temp=split("\t",$l); 
 	my $chr=$temp[0]; 
	$chr=~s/chr//g; 
	my $p1=$temp[1]+1; 
	my $p2=$temp[2]+1; 
	my $p3=$p1+1; 
	my $p4=$p2-1; 
	my @inf=split(":",$temp[3]); 
	my $intron=$inf[3]; 
	if($intron=~/^i/) {
	 my $site1=$chr.":".$p1; 
	 my $site2=$chr.":".$p2; 	
	 my $site3=$chr.":".$p3;
	 my $site4=$chr.":".$p4;
 
	 $splice_s{$site1}=1; 
	 $splice_s{$site2}=1; 
	 $splice_s{$site3}=1; 
     $splice_s{$site4}=1;
 
		}						
}

foreach my $l (`cat $f_in`) 
	{
		my $ltr=$l; 
		chomp($ltr); 
		if($ltr=~/^Hugo/ || $ltr=~/version 2\.4/) { print OUT $ltr,"\n"; }
		else { 
		my @temp=split("\t",$ltr);
	#	print $temp[8],"\n";  
		if($temp[8] eq "Intron") 
		 {
			my $chr=$temp[4]; 
			my $type=$temp[9]; 
			$chr=~s/chr//g; 
			my $ref=$temp[10]; 
			my $var=$temp[12]; 
			my $pos=$temp[5]; 
			if($type eq "SNP") 
			{
				my $p1=$chr.":".$pos; 
				if(defined $splice_s{$p1})
				{
				$ltr=~s/Intron/Splice_Site/g; 				
				} 
				print OUT $ltr,"\n"; 
			}				
			elsif($type eq "DEL")
            {
				for(my $i=0;$i<length($ref);$i++)	
                {
				my $posn=$pos+$i; 
				my $p1=$chr.":".$posn;
                if(defined $splice_s{$p1})
                {
                $ltr=~s/Intron/Splice_Site/g;  last;         
                }
				} 
                print OUT $ltr,"\n";
			}
		
            elsif($type eq "INS")
            {
				my $posn=$pos+1; 
                my $p1=$chr.":".$pos;
				my $p2=$chr.":".$posn;			
                if(defined $splice_s{$p1} && $splice_s{$p2})
                {
                $ltr=~s/Intron/Splice_Site/g;
                }
                print OUT $ltr,"\n";
            }		
		 }
		 else { print OUT $ltr,"\n"; }	
		}
	}

close IN; 
close OUT; 
