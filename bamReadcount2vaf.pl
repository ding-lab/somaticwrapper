=head1
    
    Hua Sun
    
    10/2/2019 add total depth in output
    1/2/2019
    
    perl bamReadcount_vaf.pl -s sample -l loci.txt input > output
    
    //loci.txt
    #chr start end ref alt
    
    # [ NOTE: only use the first 5 columns ]
    
    
    //input bamReadcount result
    e.g. (Denali/Katmai)
    /diskmnt/Software/bin/bam-readcount-0.7.4/mybuild/bin/bam-readcount -q 10 -b 20 -f Homo_sapiens_assembly19.fasta -l loci.txt sample.bam > bamreadcount.output
    1       898870  C       82      =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:49:58.90:32.76:36.84:6:43:0.50:0.00:4.61:7:0.14:75.86:0.49    G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  T:33:60.00:30.33:37.00:2:31:0.54:0.01:30.88:3:0.22:76.00:0.51   N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
    
    //output
    #sample chr start end ref alt   total_readcount readcound_ref   readcount_alt   vaf
    
=cut

use strict;
use Getopt::Long;

my $sample;
my $loci;
GetOptions(
    "l:s" => \$loci,
    "s:s" => \$sample
);

# from bamreadcount result
my $file = shift;
my @bamRC=`cat $file`;

# from loci
my @dataLoci=`cut -f 1-5 $loci`;

# make hash for bamReadcount data
my %hashBamRC;
my $key;
foreach (@bamRC){
    my @arr = split("\t");
    $key = "$arr[0]\t$arr[1]"; #chr,pos
    $hashBamRC{$key} = $_;
}

# extract vaf
my ($chr,$start,$end,$ref,$alt);
my ($val,$str);
foreach (@dataLoci){
   chomp;
   my ($chr,$start,$end,$ref,$alt) = split("\t");
   # print $chr,"\t",$start,"\t",$end,"\t",$ref,"\t",$alt,"\n"; 
    #<STDIN>; 
    $key = "$chr\t$start";

    if (exists $hashBamRC{$key}){
        $val = $hashBamRC{$key};
        $str = &ExtractVAF($ref,$alt,$val);
        print "$sample\t$_\t$str\n";    # output
    }
    
}


exit;

####################################################

sub ExtractVAF
{
    my ($ref,$alt,$val) = @_;

        
    my @arr = split("\t", $val);
    my $rc_ref = $arr[2];
    
    my ($r_ref, $r_alt);
    my $ref_nucleotide;
    my $alt_nucleotide;
    
    $r_ref=0;
    $r_alt=0;
    
    if ($ref eq '-'){
        $r_ref = $1 if ($val=~/\t$rc_ref\:(\d+)\:/);
        $r_alt = $1 if ($val=~/\t\+$alt\:(\d+)\:/);
    } elsif ($alt eq '-'){
        $ref_nucleotide = substr($ref, 0, 1);
        $r_ref = $1 if ($val=~/\t$ref_nucleotide\:(\d+)\:/);
        $r_alt = $1 if ($val=~/\t\-$ref\:(\d+)\:/); 
    } else {
        $ref_nucleotide = substr($ref, 0, 1);
        $alt_nucleotide = substr($alt, 0, 1);
        $r_ref = $1 if ($val=~/\t$ref_nucleotide\:(\d+)\:/);
        $r_alt = $1 if ($val=~/\t$alt_nucleotide\:(\d+)\:/);
    }
    
    # total read count
    my $total_depth = $r_ref + $r_alt;
    
    my $vaf;
    if ($total_depth==0){
        $vaf = 0;
    } else {
        $vaf = sprintf("%.2f", $r_alt/($total_depth)*100);
    }
    
    return "$total_depth\t$r_ref\t$r_alt\t$vaf";
}


