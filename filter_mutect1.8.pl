#!/usr/bin/perl

## normal <=2% 
### mutect1.8 filtering ###

## turn off minvaf cut-off, Apr 20, 2020 ##

use strict;
use warnings;
die unless @ARGV == 6;
### samtools ##

my ($samtools,$f_m,$f_filter_out,$mincov_t,$mincov_n,$minvaf)=@ARGV;

#my $f_m=$run_dir."/mutect/mutect.raw.vcf";
#my $f_filter_out=$run_dir."/mutect/mutect.filtered.vcf";

### minimum vaf for tumor 0.05, column 10 ###
## maximum vaf for normal 0.02, column 11 ###

## turn off minvaf cut-off ##

### get the bam path ##

my @temp=split(/\//,$f_m); 

my $path_d="/"; 

for(my $i=1;$i<scalar @temp-2; $i++)
{
   $path_d.=$temp[$i]."/";
}

my $f_bam_n=$path_d.$temp[-3].".N.bam";
my $f_bam_t=$path_d.$temp[-3].".T.bam";

## read absolute path ##
#my $f_bam_n_abs=`readlink -f $f_bam_n`; 
#my $f_bam_t_abs=`readlink -f $f_bam_t`; 
#chomp($f_bam_n_abs); 
#chomp($f_bam_t_abs); 
#print $f_bam_n_abs,"\n";
#print $f_bam_t_abs,"\n";
#my $last_bam_n_abs=(split(/\//,$f_bam_n_abs))[-1];
#my $last_bam_t_abs=(split(/\//,$f_bam_t_abs))[-1];
my $sn_n; 
my $sn_t; 
#foreach my $l (`$samtools view -H $f_bam_n`) 
#	{
#		my $ltr=$l; 
#		chomp($ltr); 

#		if($ltr=~/^\@RG/) { 
#			my @temp2=split("\t",$ltr); 
#			$sn_n=(split(/\:/,$temp2[-1]))[1];
#	    	last; 
#		  }
#	}


#foreach my $l (`$samtools view -H $f_bam_t`)
 #   {
  #      my $ltr=$l;
  #      chomp($ltr);
  #      if($ltr=~/^\@RG/) {
  #          my @temp2=split("\t",$ltr);
  #          $sn_t=(split(/\:/,$temp2[-1]))[1];
#			last; 
#          }
#    }

foreach my $l (`$samtools view -H $f_bam_n`)
    {
        my $ltr=$l;
        chomp($ltr);
        if($ltr=~/^\@RG/) {
            my @temp2=split("\t",$ltr);
            for(my $i=0;$i<scalar @temp2;$i++)
            {
            if($temp2[$i]=~/^SM/)
            {
            $sn_n=(split(/\:/,$temp2[$i]))[1];
            last;
            }
          }
    }
 }

foreach my $l (`$samtools view -H $f_bam_t`)
    {
        my $ltr=$l;
        chomp($ltr);
        if($ltr=~/^\@RG/) {
            my @temp2=split("\t",$ltr);
            for(my $i=0;$i<scalar @temp2;$i++)
            {
            if($temp2[$i]=~/^SM/)
            {
            $sn_t=(split(/\:/,$temp2[$i]))[1];
            last;
            }
            }
          }
    }

#print $sn_n,"\t",$sn_t,"\n";

#<STDIN>;
my $min_vaf_somatic=$minvaf;
my $max_vaf_germline=0.02;
#my $min_coverage=$mincov;
my $tumor_normal_order=-1; 

open(IN,"<$f_m");
open(OUT,">$f_filter_out");

while(<IN>) 
	{
		my $l=$_; 
		chomp($l); 
		if($l=~/^#/) { #print OUT $l,"\n";
   		if($l=~/^#CHROM/) {  
		my @temphead=split("\t",$l); 
		print OUT $temphead[0]; 
		for(my $i=1;$i<=8;$i++) 
		{
		print OUT "\t",$temphead[$i]; 
		}

		if($temphead[9] eq $sn_t && $temphead[10] eq $sn_n) { print OUT "\t","TUMOR","\t","NORMAL","\n"; $tumor_normal_order=1; }
		if($temphead[9] eq $sn_n && $temphead[10] eq $sn_t) { print OUT "\t","NORMAL","\t","TUMOR","\n"; $tumor_normal_order=0; }
		}
		else { print OUT $l,"\n"; } 
 		}

		else { 
		#print $tumor_normal_order,"\n"; 
		#<STDIN>;
		my @temp=split("\t",$l); 
		if($tumor_normal_order==-1) { last; }	
		my $tumor=$temp[9]; 
		my $normal=$temp[10];

		if($tumor_normal_order==0) { 
			$tumor=$temp[10];
			$normal=$temp[9]; 
		}

		my @tempt=split(":",$tumor); 
		my @tempn=split(":",$normal); 
		my @readt=split(",",$tempt[1]); 
		my @readn=split(",",$tempn[1]);

		#print $l,"\n";
		#print $readn[0], "\t", $readn[1],"\n";
		#<STDIN>;
		#print $readt[0], "\t", $readt[1],"\n";

#### readn[0]: read count for ref allele in normal ##
### readn[1]: read count for alt allele in normal ##
#### readt[0]: read count for ref allele in tumor ##
### readt[1]: read count for alt allele in tumor ##
 		
		my $tot_n=$readn[0]+$readn[1]; 
		my $tot_t=$readt[0]+$readt[1];										

		if($tot_n==0 || $tot_t==0) { next; }

		my $vaf_t=$readt[1]/$tot_t;
		my $vaf_n=$readn[1]/$tot_n;
		#print $vaf_t,"\t",$vaf_n,"\n";
## apply filtering ##
		if($temp[6] eq "PASS" && $vaf_n<=$max_vaf_germline && $tot_n>=$mincov_n && $tot_t>=$mincov_t) 
			{
				print OUT $l,"\n"; 
			}	
		}	
	
}


