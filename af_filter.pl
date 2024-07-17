#!/usr/bin/perl
use strict;
use warnings;

## filtering based on population af ##

my($f_vep_in,$f_vcf_in,$f_vep_out,$f_vcf_out)=@ARGV;

#my $indel_max_size=100;
my $inforheader="Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|ExAC_AF|ExAC_Adj_AF|ExAC_AFR_AF|ExAC_AMR_AF|ExAC_EAS_AF|ExAC_FIN_AF|ExAC_NFE_AF|ExAC_OTH_AF|ExAC_SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS";

open(INvep,"<$f_vep_in"); 
open(INvcf,"<$f_vcf_in");
open(OUTvep,">$f_vep_out"); 
open(OUTvcf,">$f_vcf_out");

my $gnomad_af_ind;
my %pass_var=();
my $var; 
my $info; 
my @t2;
my $gnomad_af; 

my @t=split(/\|/,$inforheader);

for(my $i=0;$i<scalar @t;$i++)
{ 
 if($t[$i] eq "gnomAD_AF")
 {
  $gnomad_af_ind=$i;   
 }
}

while(<INvep>)
	{
	  my $ltr=$_; 
      chomp($ltr);

      if($ltr=~/^#/) { 
        print OUTvep $ltr,"\n";
		}
	  else {
	    @t=split("\t",$ltr);
		$info=$t[7]; 
        @t2=split(/\|/,$info);
        if(scalar @t2>$gnomad_af_ind)
        {
        $gnomad_af=$t2[$gnomad_af_ind];
        }
        else { $gnomad_af=""; }

       # print $ltr,"\n"; 
       # print $gnomad_af,"\n";
       # <STDIN>;
        print $gnomad_af,"\n";
        if(($gnomad_af=~ /^-?\d+(\.\d+)?$/) && $gnomad_af>0.005)
        {
            next; 
        }
        else 
        {
         $var=$t[0]."_".$t[1]."_".$t[3]."_".$t[4];
         $pass_var{$var}=1; 
         print OUTvep $ltr,"\n";  
        }
		} 
	}

while(<INvcf>)
	{
	  my $ltr=$_; 
      chomp($ltr);

      if($ltr=~/^#/) { 
        print OUTvcf $ltr,"\n";
		}
	  else {
	    @t=split("\t",$ltr);
         $var=$t[0]."_".$t[1]."_".$t[3]."_".$t[4];
         if(defined $pass_var{$var}) { print OUTvcf $ltr,"\n";  }  
        } 
	}


close INvep;
close INvcf;
close OUTvep;
close OUTvcf;

