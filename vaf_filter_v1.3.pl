#!/usr/bin/perl

## tumor >= 5% and normal <=1% 

### pindel tumor >=10% since the vaf calculation underestimate the ref coverage ##

### add the filtering for indel length ##

use strict;
use warnings;
die unless @ARGV == 5;
my ($run_dir,$min_vaf_somatic,$min_coverage_t,$min_coverage_n,$indel_max_size)=@ARGV; 

my $f_m=$run_dir."/merged.withmutect.vcf"; 
my $f_filter_out=$run_dir."/merged.filtered.withmutect.vcf";
my $f_vaf_out=$run_dir."/merged.withmutect.vaf";
#my $min_vaf_somatic=0.05;
#my $min_vaf_pindel=0.1;
my $max_vaf_germline=0.02; 
#my $min_coverage=20; 
#my $indel_max_size=100; 

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

         if(length($ref)>$indel_max_size || length($var)>$indel_max_size)  { next; }
 
		 #if($info=~/strelka-varscan/) and strelka-mutect ##
		 #if($info=~/set\=sindel-varindel/ || $info=~/set\=sindel-pindel/ ##
       	if($info=~/set\=sindel-varindel/ || $info=~/set\=sindel-pindel/) 
#		if($info=~/set\=strelka-varscan/ || $info=~/set\=strelka-mutect/)
		 {

			#print $info,"\n"; 
			#<STDIN>;
			#print $ltr,"\n";
		
			$vaf_n=$temp[9];
		 	$vaf_t=$temp[10];
		 	$ref=$temp[3]; 
		 	$var=$temp[4];
		 	$r_tot=0;
			#print $ref,"\t",$var,"\n"; 
			#print $vaf_n,"\n"; 
			#print $vaf_t,"\n";

			#<STDIN>;	
		 	@temp2=split(":",$vaf_n);
			#%rc=();			
			#print $vaf_n,"\n";
			#print $temp2[0],"\t",$temp2[1],"\t",$temp2[4],"\t",$temp2[7],"\n";
			#<STDIN>;
			#my @temp3=split(",",$temp2[0]);
			#$rc{'A'}=$temp3[0];

#Somatic indel indel allele frequency is: alt_t1count / (ref_t1count + alt_t1count)
#...where:
#ref_counts = value of FORMAT column value “TAR”
#alt_counts = value of FORMAT column value “TIR”
#ref_t1count = $ref_counts[0] (use the tier1 counts -- the first value in the comma-delimited list)
#alt_t1count = $alt_counts[0] (...likewise..)

			my $rcvar=(split(",",$temp2[-2]))[0];
	 	 	my $rctot=(split(",",$temp2[-3]))[0]+(split(",",$temp2[-2]))[0]; 
			my $rcref=$rctot-$rcvar; 

		 	#$rc{'C'}=(split(",",$temp2[-2]))[0];
			#print $vaf_n,"\n";
			#print $vaf_t,"\n";	
			#<STDIN>;
#			foreach my $nt (keys %rc) 
#			{
#				$r_tot+=$rc{$nt}; 
				#print $rc{$nt},"\n";
#			}

			@temp2=split(":",$vaf_t);

#Somatic indel indel allele frequency is: alt_t1count / (ref_t1count + alt_t1count)
#...where:
#ref_counts = value of FORMAT column value “TAR”
#alt_counts = value of FORMAT column value “TIR”
#ref_t1count = $ref_counts[0] (use the tier1 counts -- the first value in the comma-delimited list)
#alt_t1count = $alt_counts[0] (...likewise..)

	        my $rc2var=(split(",",$temp2[-2]))[0];
            my $rc2tot=(split(",",$temp2[-3]))[0]+(split(",",$temp2[-2]))[0];
			my $rc2ref=$rc2tot-$rc2var; 		

			#print $rcref,"\t",$rcvar,"\t",$rctot,"\n";
			#print $rc2ref,"\t",$rc2var,"\t",$rc2tot,"\n";
			#<STDIN>;

		#print $rc{$ref},"\t",$rcvar,"\t",$rc2{$ref},"\t",$rc2var,"\n"; 
		#<STDIN>;
			print OUT2 $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t",$info,"\t",$rcref,"\t",$rcref/$rctot,"\t",$rcvar,"\t",$rcvar/$rctot,"\t",$rc2ref,"\t",$rc2ref/$rc2tot,"\t",$rc2var,"\t",$rc2var/$rc2tot,"\n"; 

			if($rc2var/$rc2tot>=$min_vaf_somatic && $rcvar/$rctot<=$max_vaf_germline && $rc2tot>=$min_coverage_t && $rctot>=$min_coverage_n) 
			{
				$ltr=~s/SVTYPE=//g;
		        print OUT1 $ltr,"\n";
			} 	
		}


#Somatic indel indel allele frequency is: alt_t1count / (ref_t1count + alt_t1count) for strelka

#...where:

#ref_counts = value of FORMAT column value “TAR”
#alt_counts = value of FORMAT column value “TIR”
#ref_t1count = $ref_counts[0] (use the tier1 counts -- the first value in the comma-delimited list)
#alt_t1count = $alt_counts[0] (...likewise..)
         if($info=~/set\=strelka-varscan/ || $info=~/set\=strelka-mutect/)
       #if($info=~/set\=sindel-varindel/ || $info=~/set\=sindel-pindel/)
          	{

            #print $info,"\n"; 
            #<STDIN>;
#            print $ltr,"\n";

            $vaf_n=$temp[9];
            $vaf_t=$temp[10];
            $ref=$temp[3];
            $var=$temp[4];
            $r_tot=0;
 #           print $ref,"\t",$var,"\n";
  #          print $vaf_n,"\n";
   #         <STDIN>;
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
            #$rc{'A'}=(split(",",$temp2[2]))[0];
            #$rc{'C'}=(split(",",$temp2[3]))[0];
            #$rc{'G'}=(split(",",$temp2[6]))[0];
            #$rc{'T'}=(split(",",$temp2[9]))[0];
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

            #$rc2{'A'}=(split(",",$temp2[2]))[0];
            #$rc2{'C'}=(split(",",$temp2[3]))[0];
            #$rc2{'G'}=(split(",",$temp2[6]))[0];
            #$rc2{'T'}=(split(",",$temp2[9]))[0];

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

        	if($rc2var/$r_tot2>=$min_vaf_somatic && $rcvar/$r_tot<=$max_vaf_germline && $r_tot2>=$min_coverage_t && $r_tot>=$min_coverage_n)
            {
                print OUT1 $ltr,"\n";
            }
        }



		# called snv called by varscan and mutect, indel called varscan and strelka, 	
		if($info=~/set\=varscan-mutect/ || $info=~/set\=varindel-sindel/ || $info=~/set\=varindel-pindel/)
		{
		   	$vaf_n=$temp[9];
        	$vaf_t=$temp[10];
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
		if($tdp_var/($tdp_var+$tdp_ref) >=$min_vaf_somatic && $ndp_var/($ndp_var+$ndp_ref)<=$max_vaf_germline && $tdp_var+$tdp_ref>=$min_coverage_t && $ndp_var+$ndp_ref>=$min_coverage_n) 	
			{
			$ltr=~s/SVTYPE=//g;
			print OUT1 $ltr,"\n"; 
			}
		}



	}
}
