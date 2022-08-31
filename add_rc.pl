### add read counts and vafs for ref and var to maf file ##
### Song Cao ##

#!/usr/bin/perl
use strict;
use warnings;
(my $usage = <<OUT) =~ s/\t+//g;
This script will add readcounts for ref and var to maf file 
perl run_dir f_maf f_out
OUT

die $usage unless @ARGV == 3;
my ($run_dir,$f_maf,$f_out)=@ARGV;
my $sn_1=""; 
my $sn_2=""; 
my %vaf_rc=(); 

open(OUT,">$f_out");

foreach my $l (`cat $f_maf`) 
	{

	my $ltr=$l; 
	chomp($ltr);
	my @temp=split("\t",$ltr); 
 
	if($ltr=~/^Hugo/) { print OUT $ltr,"\n"; }

	else {

	 	$sn_1=$temp[15]; $sn_1=~s/_T$//g; 

		if($sn_2 eq "" || $sn_2 ne $sn_1)
		{

		 %vaf_rc=();  
		 $sn_2=$sn_1; 

		 my $f_vaf=$run_dir."/".$sn_1."/merged.withmutect.vaf"; 
		 open(IN,"<$f_vaf"); 
		 my $id; 

		 my $n_ref;
		 my $n_var; 
		 my $t_ref; 
		 my $t_var; 

		 while(<IN>)
		  {
		    my $l2=$_; 
			chomp($l2); 
			my @temp2=split("\t",$l2); 
			my $ref=$temp2[3]; 
			my $var=$temp2[4];
			my $equal=0;
			#print $ref,"\t",$var,"\n";
			for(my $i=0;$i<length($ref) && $i<length($var); $i++)
			{
			 if(substr($ref,$i,1) eq substr($var,$i,1)) 
			 {
				$equal++;
			 }
			} 
			# insertion ##
			my $shift=$equal; 
			if($equal eq length($ref)) { $shift=0; }

			my $pos1=$temp2[1]+$shift; 
		
		    my $id=$temp2[0]."_".$pos1;
 
		  	$n_ref=$temp2[6];
			$n_var=$temp2[8];
			$t_ref=$temp2[10]; 
			$t_var=$temp2[12]; 
			#print $sn_1,"\t",$n_ref,"\t",$n_var,"\t",$t_ref,"\t",$t_var,"\n"; 
			#<STDIN>;
        	$vaf_rc{$id}=$n_ref."_".$n_var."_".$t_ref."_".$t_var; 	  	 	
		  }
		  
         close IN; 
		}	

		# add read counts for tumor and normal ##	    
	## t_depth	t_ref_count	t_alt_count	n_depth	n_ref_count	n_alt_count
	## 40-45
	   my $pos2=$temp[5];

	  # if($temp[12] eq "-")
	#	{
	#		$pos-=1;
	#	}
	
		my $id2=$temp[4]."_".$pos2;
		if(defined $vaf_rc{$id2}) 
		{ 
		my $vaf=$vaf_rc{$id2};

		#print $id2,"\n";
	#	print $vaf,"\n";
	#	print $sn_1,"\n";
		#<STDIN>; 
		my @temp3=split("_",$vaf); 
	   	my $t_depth=$temp3[2]+$temp3[3]; 
		my $t_ref_count=$temp3[2]; 
		my $t_alt_count=$temp3[3];
		  								
	    my $n_depth=$temp3[0]+$temp3[1];
        my $n_ref_count=$temp3[0];
        my $n_alt_count=$temp3[1];
	
		print OUT $temp[0]; 

		for(my $i=1;$i<39;$i++)
		{ 
		print OUT "\t",$temp[$i]; 
		}

		print OUT "\t",$t_depth,"\t",$t_ref_count,"\t",$t_alt_count,"\t",$n_depth,"\t",$n_ref_count,"\t",$n_alt_count; 						
		for(my $i=45;$i<scalar @temp;$i++) 
		{
		print OUT "\t",$temp[$i]; 		
		}

		print OUT "\n"; 
		
		}	
		
		else { 

		  print $id2,"\n";
		  print $sn_1,"\n";
			
			}
	}
	}
