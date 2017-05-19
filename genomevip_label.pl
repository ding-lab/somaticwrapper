#!/usr/bin/env perl
# ---------------------
# @name GenomeVIP utility to annotate with variant calling tool program
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# @version
#
# @syntax  anno_genomevip.pl  program_name   source.vcf  destination.vcf
# ---------------------
use strict;
use warnings;

my ($progid, $infile, $outfile, @header, $encoding, %assoc);


$progid = $ARGV[0]; $infile = $ARGV[1];  $outfile = $ARGV[2];
%assoc = ("VarScan"=>1,"Pindel"=>2,"BreakDancer"=>3,"GenomeSTRiP"=>4,"Strelka"=>5,"GATK"=>6,"MuTect"=>7);
$encoding="VarScan,1;Pindel,2;BreakDancer,3;GenomeSTRiP,4;Strelka,5;GATK,6;MuTect,7";

# Insert custom info line after last existing info line, then update info
open (IN, "<", $infile) || die "Error: cannot open input file $infile"; close(IN);
@header = `grep '^#' $infile`;
my $idx=-1;
for (my $i=scalar(@header)-1; $i>= 0; $i--) {
    if ($header[$i] =~ /^##INFO=/) {
	$idx = $i;
	last;
    }
}
open (OUT,">", $outfile) || die "Error: cannot open output file $outfile";
for (my $i=0; $i<scalar(@header); $i++) {
    print OUT $header[$i];
    print OUT  "##INFO=<ID=GENOMEVIP,Number=.,Type=Integer,Description=\"Call support from programs ".$encoding."\">\n" if ($i==$idx);
}
open (IN, "<", $infile);
while(<IN>) {
    if($_ !~ /^#/) {
	chomp; 
	my @a = split /\t/;
	$a[7] .= ";GENOMEVIP=". $assoc{$progid};
	print OUT join("\t", @a)."\n";
    }
}
close(IN);
close(OUT);
