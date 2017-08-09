#!/usr/bin/perl

use strict;
use warnings;

my @dir=("/gscuser/scao/gc2521/dinglab/cptac_prospective_samples/exome/somatic/CO", "/gscuser/scao/gc2521/dinglab/cptac_prospective_samples/exome/somatic/CO2", "/gscuser/scao/gc2521/dinglab/cptac_prospective_samples/exome/somatic/BRCA2", "/gscuser/scao/gc2521/dinglab/cptac_prospective_samples/exome/somatic/BRCA3", "/gscuser/scao/gc2521/dinglab/cptac_prospective_samples/exome/somatic/OV");
#my @dir=("/gscuser/scao/gc2521/dinglab/cptac_prospective_samples/exome/somatic/CO2");

my $tool_dir="/gscmnt/gc2741/ding/qgao/tools/vcf2maf-1.6.11";
foreach my $d (@dir)
{
#	print "$d\n";
	foreach my $s (glob("$d/*"))
	{
		my @samp = split(/\//,$s);
		system("ln -s $s/merged.vcf $samp[-1].vcf");
		system("ln -s $s/merged.VEP.vcf $samp[-1].vep.vcf");
		system("bsub -oo $samp[-1].log perl $tool_dir/vcf2maf.pl --input-vcf $samp[-1].vcf --output-maf $samp[-1].maf --tumor-id $samp[-1]\_T --normal-id $samp[-1]\_N --ref-fasta /gscuser/scao/gc3027/fasta/GRCh37V1/GRCh37-lite-chr_with_chrM.fa --filter-vcf $tool_dir/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz");
	}
}

