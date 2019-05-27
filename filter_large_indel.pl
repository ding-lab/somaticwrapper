#!/usr/bin/perl
use strict;
use warnings;

## filter large indel (>100) in pindel output ##
## remove END pos in pindel vcf output, which can cause gatk merging problems ##

my ($f_indel_in,$f_indel_out,$indel_max_size)=@ARGV;

#my $indel_max_size=100;

open(IN,"<$f_indel_in"); 
open(OUT,">$f_indel_out"); 

while(<IN>)
	{
	  my $ltr=$_; 
      chomp($ltr);
      if($ltr=~/^#/) { 
        if($ltr=~/^#CHROM/) {
        my @temphead=split("\t",$ltr);
        print OUT $temphead[0];
		
        for(my $i=1;$i<=8;$i++)
        {
          print OUT "\t",$temphead[$i];
        }

        if($temphead[9]=~/\.T$/ && $temphead[10]=~/\.N$/) { print OUT "\t","TUMOR","\t","NORMAL","\n"; last; }
        if($temphead[9]=~/\.N$/ && $temphead[10]=~/\.T$/) { print OUT "\t","NORMAL","\t","TUMOR","\n";  }

        }
        else { print OUT $ltr,"\n"; }
		}
	  else {
	    my @temp=split("\t",$ltr);
		my $ref=$temp[3];
        my $var=$temp[4];
		my $info=$temp[7]; 
		$info=~s/END=\d+;//g;  
        if(length($ref)>$indel_max_size || length($var)>$indel_max_size)  { next; }
	    else
		{
			print OUT $temp[0];
			for(my $i=1;$i<=6;$i++) 
			{
				print OUT "\t",$temp[$i];   
			} 
			print OUT "\t",$info;
		    for(my $i=8;$i<=10; $i++)  
            {
                print OUT "\t",$temp[$i];   
            } 
			print OUT "\n"; 
		}
		} 
	}

close OUT;
close IN; 
#my $f_m=$run_dir."/mutect/mutect.raw.vcf";
#my $f_filter_out=$run_dir."/mutect/mutect.filtered.vcf";

### minimum vaf for tumor 0.05, column 10 ###
## maximum vaf for normal 0.02, column 11 ###
## minimum coverage 20 ###

### get the bam path ##
