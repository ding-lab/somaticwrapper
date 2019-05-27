#!/usr/bin/perl

## tumor >= 5% and normal <=1% 

### pindel tumor >=10% since the vaf calculation underestimate the ref coverage ##

### add the filtering for indel length ##

use strict;
use warnings;
die unless @ARGV == 1;
my ($run_dir)=@ARGV; 

my $f_m=$run_dir."/merged.withmutect.vcf"; 
my $f_filter_out=$run_dir."/merged.filtered.withmutect.vcf";
my $f_vaf_out=$run_dir."/merged.withmutect.vaf";
my $min_vaf_somatic=0.05;
my $min_vaf_pindel=0.1;
my $max_vaf_germline=0.02; 
my $min_coverage=20; 
my $indel_max_size=100; 

open(OUT1,">$f_filter_out");
open(OUT2,">$f_vaf_out"); 

foreach my $l (`cat $f_m`) 
	{
		my $ltr=$l;
		chomp($ltr);  
		if($ltr=~/^#/) { print OUT1 $ltr,"\n"; next;  }
		else {
		 my @temp=split("\t",$ltr); 
		 #my $ref=$temp[3];
		 #my $var=$temp[4];


	 	 my $info=$temp[7];
		 my @temp2; 	
		 my %rc; 
		 my %rc2;
		 my $r_tot;  
		 my $r_tot2; 
		 my $vaf_n;
  		 my $vaf_t;
		 my $ref;
		 my $var;
		 my $nt;
		 my $ndp_ref;
		 my $ndp_var;
		 my $tdp_ref;
		 my $tdp_var; 

         $ref=$temp[3];
         $var=$temp[4];

         if(length($ref)>=$indel_max_size || length($var)>=$indel_max_size)  { next; }
 
		 #if($info=~/strelka-varscan/) 
		if($info=~/set\=strelka-varscan/ || $info=~/set\=strelka-mutect/ || $info=~/set\=pindel-sindel/)
		 {
			#print $info,"\n"; 
			#<STDIN>;
	
			$vaf_n=$temp[11];
		 	$vaf_t=$temp[12];
		 	$ref=$temp[3]; 
		 	$var=$temp[4];
		 	$r_tot=0; 
	
		 	@temp2=split(":",$vaf_n);
			%rc=();			
			#print $vaf_n,"\n";
			#print $temp2[0],"\t",$temp2[1],"\t",$temp2[4],"\t",$temp2[7],"\n";
			#<STDIN>;
			#my @temp3=split(",",$temp2[0]);
			#$rc{'A'}=$temp3[0];

	 	 	$rc{'A'}=(split(",",$temp2[0]))[0]; 
		 	$rc{'C'}=(split(",",$temp2[1]))[0];
		 	$rc{'G'}=(split(",",$temp2[4]))[0];
		 	$rc{'T'}=(split(",",$temp2[7]))[0];
			#print $vaf_n,"\n";
			#print $vaf_t,"\n";	
			#<STDIN>;
			foreach my $nt (keys %rc) 
			{
				$r_tot+=$rc{$nt}; 
				#print $rc{$nt},"\n";
			}

			@temp2=split(":",$vaf_t);

        	%rc2=();

         	$rc2{'A'}=(split(",",$temp2[0]))[0];
         	$rc2{'C'}=(split(",",$temp2[1]))[0];
         	$rc2{'G'}=(split(",",$temp2[4]))[0];
         	$rc2{'T'}=(split(",",$temp2[7]))[0];

        	foreach $nt (sort keys %rc2)
        	{
            	$r_tot2+=$rc2{$nt};
				#print $rc2{$nt},"\n";
        	}
		#print $ltr,"\n";
		my @vars=split(",",$var); 
		my $rcvar=0;
		my $rc2var=0; 
		foreach my $v (@vars)
		{
		 $rcvar+=$rc{$v};
		 $rc2var+=$rc2{$v}; 
		}		
		
		#print $rc{$ref},"\t",$rcvar,"\t",$rc2{$ref},"\t",$rc2var,"\n"; 
		#<STDIN>;
		print OUT2 $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t",$info,"\t",$rc{$ref},"\t",$rc{$ref}/$r_tot,"\t",$rcvar,"\t",$rcvar/$r_tot,"\t",$rc2{$ref},"\t",$rc2{$ref}/$r_tot2,"\t",$rc2var,"\t",$rc2var/$r_tot2,"\n"; 

		if($rc2var/$r_tot2>=$min_vaf_somatic && $rcvar/$r_tot<=$max_vaf_germline && $r_tot2>=$min_coverage && $r_tot>=$min_coverage) 
			{
		        print OUT1 $ltr,"\n";
			} 	
		}
	
		elsif($info=~/set\=varscan-mutect/ || $info=~/set\=varindel-sindel/ || $info=~/set\=pindel-varindel/ )
		{
		   	$vaf_n=$temp[11];
        	$vaf_t=$temp[12];
			@temp2=split(":",$vaf_n); 
			#print $vaf_n,"\n";
			my @ndp4=split(",",$temp2[3]);

			if(scalar @ndp4<4) { @ndp4=split(",",$temp2[2]);  }  
			$ndp_ref=$ndp4[0]+$ndp4[1];
			$ndp_var=$ndp4[2]+$ndp4[3];
			#print $vaf_t,"\n";
            @temp2=split(":",$vaf_t);
            my @tdp4=split(",",$temp2[3]);
            if(scalar @tdp4<4) { @tdp4=split(",",$temp2[2]);  }	

            $tdp_ref=$tdp4[0]+$tdp4[1];
            $tdp_var=$tdp4[2]+$tdp4[3];
	        print OUT2 $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t",$info,"\t",$ndp_ref,"\t",$ndp_ref/($ndp_ref+$ndp_var),"\t",$ndp_var,"\t",$ndp_var/($ndp_var+$ndp_ref),"\t",$tdp_ref,"\t",$tdp_ref/($tdp_ref+$tdp_var),"\t",$tdp_var,"\t",$tdp_var/($tdp_var+$tdp_ref),"\n";  
		if($tdp_var/($tdp_var+$tdp_ref) >=$min_vaf_somatic && $ndp_var/($ndp_var+$ndp_ref)<=$max_vaf_germline && $tdp_var+$tdp_ref>=$min_coverage && $ndp_var+$ndp_ref>=$min_coverage) 	
			{
			$ltr=~s/SVTYPE=//g;
			print OUT1 $ltr,"\n"; 
			}
		}

        elsif($info=~/set\=pindel/)
          {

			$vaf_t=$temp[10];
            $vaf_n=$temp[9];
			if(!($vaf_t=~/\:/)) { next; } 
			if(!($vaf_n=~/\:/)) { next; } 

			@temp2=split(":",$vaf_n);	
			my @ndp2=split(",",$temp2[1]);
            $ndp_ref=$ndp2[0];
            $ndp_var=$ndp2[1];
            @temp2=split(":",$vaf_t);
            my @tdp2=split(",",$temp2[1]);
            $tdp_ref=$tdp2[0];
            $tdp_var=$tdp2[1];
			#print $ndp_ref,"\t",$ndp_var,"\t",$tdp_ref,"\t",$tdp_var,"\n";
			#<STDIN>;

			print OUT2 $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t",$info,"\t",$ndp_ref,"\t",$ndp_ref/($ndp_ref+$ndp_var),"\t",$ndp_var,"\t",$ndp_var/($ndp_var+$ndp_ref),"\t",$tdp_ref,"\t",$tdp_ref/($tdp_ref+$tdp_var),"\t",$tdp_var,"\t",$tdp_var/($tdp_var+$tdp_ref),"\n";

		if($tdp_var/($tdp_var+$tdp_ref)>=$min_vaf_pindel && $ndp_var/($ndp_ref+$ndp_var)<=$max_vaf_germline && $tdp_var+$tdp_ref>=$min_coverage && $ndp_var+$ndp_ref>=$min_coverage) 
		{
			#print $ltr,"\n";    
 		$ltr=~s/SVTYPE=//g;		 
  		print OUT1 $ltr,"\n";	
		}
		
	}


	}
}
