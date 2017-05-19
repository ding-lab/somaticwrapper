######### Song Cao###########
## combine somatic variant result ##
#  combine_variant.pl #

### updated date: 04/18/2017 ###

#!/usr/bin/perl
use strict;
use warnings;

die unless @ARGV == 1;
my ($run_dir)=@ARGV; 

my @temp = split("/", $dir);
my $run_name = pop @temp;
my $f_out=$run_dir."/Analysis_Summary_Somatic_$run_name.tsv";
open(OUT,">$f_out"); 

foreach my $dir_s (`ls $run_dir`) 
	{
		my $dir_s_f=$run_dir."/".$dir_s; 
		if(-d $dir_s_f) 
			{

				my %mut_inf=();
				my %snv_varscan_tag=();
				my %indel_varscan_tag=();
				my %snv_strelka_tag=();
				my %indel_strelka_tag=();

				my $f_varscan_snv=$dir_s_f."/varscan/varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf";
				my $f_varscan_indel=$dir_s_f."/varscan/varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf";	
				my $f_strelka_snv=$dir_s_f."/strelka/strelka_out/results/strelka.somatic_snv.current_final.gvip.Somatic.VEP.vcf";
				my $f_strelka_indel=$dir_s_f."/strelka/strelka_out/results/strelka.somatic_indel.current_final.gvip.Somatic.VEP.vcf";
				foreach my $v (`cat $f_varscan_snv`) 
				{
					my $vtr=$v; 
					chomp($vtr); 
					if($vtr=~/^#/)
					{
						next; 	
					}	
					my @temp=split("\t",$vtr); 
					my $mut=$temp[0]."_".$temp[1]."_".$temp[3]."_".$temp[4]; 
					$mut_inf{$mut}=$vtr; 
															
				}					 
			}
	}
