### add read counts and vafs for ref and var to maf file ##
### Song Cao ##

#!/usr/bin/perl
use strict;
use warnings;
(my $usage = <<OUT) =~ s/\t+//g;
Get caller information
perl run_dir f_maf f_out
OUT

die $usage unless @ARGV == 3;

my ($run_dir,$f_maf,$f_out)=@ARGV;
my $sn_1=""; 
my $sn_2=""; 
my %caller=(); 

my @allchr=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"); 
my %chrlist=();  
foreach my $c (@allchr) 
	{
	$chrlist{$c}=1; 
		
	}
open(OUT,">$f_out");

foreach my $l (`cat $f_maf`) 
	{

	my $ltr=$l; 
	chomp($ltr);
	my @temp=split("\t",$ltr); 
	my $i;  
	if($ltr=~/^Hugo/) { for($i=0;$i<=44;$i++) { print OUT $temp[$i],"\t"; } print OUT "callers"; for($i=45;$i<scalar @temp;$i++) { print OUT "\t",$temp[$i]; } print OUT "\n"; }

	else {

		my $chr=$temp[4]; 
		$chr=~s/chr//g; 
		if(!defined $chrlist{$chr}) { next; } 
				
	 	$sn_1=$temp[15]; $sn_1=~s/_T//g; 

		if($sn_2 eq "" || $sn_2 ne $sn_1)
		{

		 %caller=();  
		 $sn_2=$sn_1; 
		 print $sn_1,"\n";
		 my $f_vcf=$run_dir."/".$sn_1."/merged.filtered.withmutect.vcf"; 
		 open(IN,"<$f_vcf"); 
		 my $id; 
		 #my $caller;
 
		 #my $n_ref;
		 #my $n_var; 
		 #my $t_ref; 
		 #my $t_var; 

		 while(<IN>)
		  {
		    my $l2=$_; 
			chomp($l2); 
			#print $l2,"\n";

			if($l2=~/^#/) { next; }
			#my $temp7=$temp[7];

			my @temp2=split("\t",$l2); 
			my $info2=$temp2[7];

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

			#print $info2,"\n"; <STDIN>;
			my @temp3=split(";",$info2); 
			my $info3=$temp3[-1]; 
		    $info3=~s/set=//g;
			print $id, "\t",$info3,"\n";
         	$caller{$id}=$info3; 	  	 	
		  }
 
         close IN; 
		}	

		# add read counts for tumor and normal ##	    
	## t_depth	t_ref_count	t_alt_count	n_depth	n_ref_count	n_alt_count
	## 40-44
	   my $pos2=$temp[5];

	  # if($temp[12] eq "-")
	#	{
	#		$pos-=1;
	#	}
	
		my $id2=$temp[4]."_".$pos2;
		print $id2,"\n"; 
		#<STDIN>;

		if(defined $caller{$id2}) 
		{ 

		my $callerinf=$caller{$id2};
		print $callerinf,"\n"; 
		#<STDIN>;		
		#print $id2,"\n";
	#	print $vaf,"\n";
	#	print $sn_1,"\n";
		#<STDIN>; 
		#my @temp3=split("_",$vaf); 
	   	#my $t_depth=$temp3[2]+$temp3[3]; 
		#my $t_ref_count=$temp3[2]; 
		#my $t_alt_count=$temp3[3];
		  								
	    #my $n_depth=$temp3[0]+$temp3[1];
        #my $n_ref_count=$temp3[0];
        #my $n_alt_count=$temp3[1];
	
		print OUT $temp[0]; 

		for(my $i=1;$i<=44;$i++)
		{ 
		print OUT "\t",$temp[$i]; 
		}

		print OUT "\t",$callerinf; 	

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
