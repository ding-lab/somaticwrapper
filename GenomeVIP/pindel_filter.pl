#!/usr/bin/env perl
#--------------------------------------
# @name GenomeVIP Pindel filters
# @author Beifang Niu
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# 
# @version 0.6 (maw): modularize to allow for improved debugging
# @version 0.5 (rjm): improve documentation
# @version 0.4 (rjm): (optional) user-specified log file for appending; mode not required to be specified when not performing filtering
# @version 0.3 (rjm): generalize filter hierarchy filename handling and simply code structure
# @version 0.2 (rjm): added coverage, germline, trio, minimal pool filtering; revised approach; added failed pass. Adjusted filename and parameters names; allow for commented lines
# @version 0.1 (bn):  original somatic filter, written for (tumor,normal) column order
#--------------------------------------
use strict;
use warnings;

use Cwd;
use Carp;
use FileHandle;
use IO::File;
use Getopt::Long;
use POSIX qw( WIFEXITED );
use File::Temp qw/ tempfile /;
use File::Basename;

# List of parameters accepted:
# * variants_file
# * zero_variant_support
# * apply_filter
# * mode
# * min_coverages
# * min_var_allele_freq
# * require_balanced_reads
# * remove_complex_indels
# * child_var_allele_freq
# * parents_max_num_supporting_reads
# * REF
# * date
# * heterozyg_min_var_allele_freq
# * homozyg_min_var_allele_freq
# * max_num_homopolymer_repeat_units
# * skip_filter1:           skip CvgVafStrandFilter
# * skip_pindel2vcf:        do not run pindel2vcf
# * skip_filter2:           skip Homopolymer filter 

sub CvgVafStrandFilter {
    my $infn = shift;   
    my $passfn = shift; 
    my $failfn = shift; 
    my $mode = shift;   
    my $min_coverages = shift;
    my $min_var_allele_freq = shift;
    my $require_balanced_reads = shift; # boolean
    my $remove_complex_indels = shift;  # boolean
    my $child_var_allele_freq = shift;  
    my $parents_max_num_supporting_reads = shift; 
    my $zero = shift;

# Optional filter, part 1
    my $input_fh        = IO::File->new( "$infn"         ) or die "Could not open $infn for reading $! ";
    my $filter1_pass_fh = IO::File->new( "$passfn",  ">" ) or die "Could not create $passfn for writing $!";
    my $filter1_fail_fh = IO::File->new( "$failfn",  ">" ) or die "Could not create $failfn for writing $!";

# Pindel column output
# Ref: http://seqanswers.com/forums/showthread.php?t=41121
# Ref: http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
#
# In terms of 0-based numbering for three samples:
#SampleID RefSupportingLeft RefSupportingRight AltSupportingLeft AltSupportingLeftUnique AltSupportingRight AltSupportingRightUnique
#   31          32               33                 34                    35                    36                  37
#   38          39               40                 41                    42                    43                  44
#   45          46               47                 48                    49                    50                  51


# Pindel's genotyping code takes the total reference support as max($t[32],$t[33]), etc. It may be that pindel is counting 
# lefts and rights separately so as to avoid misrepresenting total depth at a given location where left and right read fragments 
# happen to overlap at that location. In the VAF calculations below, we implicity use ref support = avg($t[32],$t[33]).

# Germline filtering options:
# minimum coverage, VAF threshold, reads are considered balanced as long as there is nonzero read support in both directions
    if ($mode eq "germline") {
        while (<$input_fh>) {
            chomp; 
            my @t = split /\s+/;
            if(  ($t[32] + $t[34] + $t[36] <  $min_coverages) || ($t[33] + $t[34] + $t[36] <  $min_coverages)  ) {
                $filter1_fail_fh->print($_."\n");
                next;
            }
            if( ($t[34] + $t[36] + $t[34] + $t[36])/($t[32] + $t[33] + $t[34] + $t[36] + $t[34] + $t[36] ) <  $min_var_allele_freq ){
                $filter1_fail_fh->print($_."\n");
                next;
            }
            if($require_balanced_reads) {  
                if ( $t[34] == 0 ||  $t[36] == 0 ) {
                    $filter1_fail_fh->print($_."\n");
                    next;
                }
            }
            $filter1_pass_fh->print($_."\n");
        }
    }

# Somatic filtering options:
# This calculation assumes sample column order is tumor/normal, as done in GenomeVIP.
# minimum coverage met for both samples; VAF threshold in tumor; zero variant
# support in normal; reads are considered balanced as long as there is nonzero
# read support in both directions in tumor
    if ($mode eq "somatic") {
        while (<$input_fh>) {
            chomp; 
            my @t = split /\s+/;
            if(  ($t[32] + $t[34] + $t[36] <  $min_coverages) || ($t[33] + $t[34] + $t[36] < $min_coverages) || ($t[39] + $t[41] + $t[43] < $min_coverages) || ($t[40] + $t[41] + $t[43] < $min_coverages)) {
                $filter1_fail_fh->print($_."\n");
                next;
            }
            if( ($t[34] + $t[36] + $t[34] + $t[36])/($t[32] + $t[33] + $t[34] + $t[36] + $t[34] + $t[36] ) > $zero ||  ($t[41] + $t[43] + $t[41] + $t[43])/($t[39] + $t[40] + $t[41] + $t[43] + $t[41] + $t[43]) < $min_var_allele_freq) {
                $filter1_fail_fh->print($_."\n");
                next;
            }
            if($require_balanced_reads) {  
                if ( $t[41] == 0 || $t[43] == 0 ) {
                    $filter1_fail_fh->print($_."\n");
                    next;
                }
            }
            if($remove_complex_indels) {  
                if ( $t[1] eq "I" || $t[1] eq "D") {
                    if ( $t[1] eq "I" || ($t[1] eq "D" && $t[4] == 0) ) {
# print "Indel filter: passed\n";
                        $filter1_pass_fh->print($_."\n");
                    } else {
                        $filter1_fail_fh->print($_."\n");
                    }
                    next;
                }
            }
            $filter1_pass_fh->print($_."\n");
        }
    }

# Trio filtering options:
# This calculation assumes sample column order is parent/parent/child, as done in GenomeVIP.
# minimum coverage met for all samples; VAF threshold in child; maximum allowed
# variant support in parents combined; reads are considered balanced as long as
# there is nonzero read support in both directions in child
    if ($mode eq "trio") {   # trio
        while (<$input_fh>) {
            chomp; 
            my @t = split /\s+/;
            if(  ($t[32] + $t[34] + $t[36] <  $min_coverages) || ($t[33] + $t[34] + $t[36] <  $min_coverages)  || ($t[39] + $t[41] + $t[43] < $min_coverages) ||  ($t[40] + $t[41] + $t[43] <  $min_coverages) ||   ($t[46] + $t[48] + $t[50] < $min_coverages) ||  ($t[47] + $t[48] + $t[50] <  $min_coverages) ) {
                $filter1_fail_fh->print($_."\n");
                next;
            }
            if( ($t[48] + $t[50] + $t[48] + $t[50])/($t[46] + $t[47] + $t[48] + $t[50] + $t[48] + $t[50]) < $child_var_allele_freq    ) { 
                $filter1_fail_fh->print($_."\n");
                next;
            }
            if( ($t[34] + $t[36]) + ($t[41] + $t[43]) >  $parents_max_num_supporting_reads  ) {
                $filter1_fail_fh->print($_."\n");
                next;
            }
            if($require_balanced_reads) {  
                if ( $t[48] == 0 ||  $t[50] == 0 ) {
                    $filter1_fail_fh->print($_."\n");
                    next;
                }
            }
            $filter1_pass_fh->print($_."\n");
        }
    }

    $filter1_fail_fh->close;
    $filter1_pass_fh->close; 
    $input_fh->close; 
}

sub runPindel2VCF {
    my $thisdir = shift;
    my $var_src = shift;
    my $var_dest = shift;
    my $ref_base = shift;
    my $logfile = shift; 
    # below are paras
    my $pindel2vcf = shift;
    my $ref = shift;
    my $date = shift;
    my $heterozyg_min_var_allele_freq = shift;
    my $homozyg_min_var_allele_freq = shift;

    # NOTE: pindel -co option seems to work on vcf output but warning messages may happend regardless
    my $pindel2vcf_command = "$pindel2vcf -r $ref -R $ref_base -p $var_src -d $date -v $var_dest -he $heterozyg_min_var_allele_freq -ho $homozyg_min_var_allele_freq  >> $thisdir/$logfile 2>&1";
    print $pindel2vcf_command."\n";
    my $return_code = system( $pindel2vcf_command );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

sub HomopolymerFilter {
    my $infn = shift;   
    my $passfn = shift; 
    my $failfn = shift; 
    my $max_num_homopolymer_repeat_units = shift;

    my $input_fh        = IO::File->new( "$infn"        ) or die "Could not open $infn for reading $! ";
    my $filter2_pass_fh = IO::File->new( "$passfn", ">" ) or die "Could not create $passfn for writing $!";
    my $filter2_fail_fh = IO::File->new( "$failfn", ">" ) or die "Could not create $failfn for writing $!";

    while ( <$input_fh> ) {
        if ( /^#/ ) { $filter2_pass_fh->print($_); next };
        my @a= split /\t/; 
        my @b = split/\;/, $a[7]; 
        for ( my $i=0; $i<scalar(@b); $i++) { 
            if ( $b[$i]=~/^HOMLEN/ ) { 
                my @c = split/=/, $b[$i]; 
                if ( $c[1] <= $max_num_homopolymer_repeat_units ) { 
                    $filter2_pass_fh->print($_); 
                }  else {
                    $filter2_fail_fh->print($_); 
                }
                last;
            } 
        }
    }
    $filter2_fail_fh->close;
    $filter2_pass_fh->close;
    $input_fh->close;
}

my $thisdir;
$thisdir=`dirname $ARGV[0]`;
chomp $thisdir;

# get paras from config file
my %paras; 
map { chomp; if ( !/^[#;]/ &&  /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
if( $paras{'apply_filter'} eq "true" &&  !exists($paras{'mode'}) ) { die "Could not detect a filtering mode for filtering !!!\n"; }
# file exist ? 
unless ( -e $paras{'variants_file'} ) { die "input indels not exist !!! \n"; }

my $zero = 0.001;
if( exists $paras{'zero_variant_support'}) {
    $zero=$paras{'zero_variant_support'};
    print("Setting zero_variant_support: $zero\n");
} else {
    print("Using default value of zero_variant_support: $zero\n");
}

my $var_file        = $paras{'variants_file'};

# Filters for coverages, vaf, and balanced reads
# Allow use case where filter 1, pindel2vcf, and filter 2 can be skipped for e.g. debugging
# Ideally all output filenames would be determined outside of this script, but for now we preserve backward
# functionality by keeping filenames unchanged but adding paramerters skip_filter1, skip_pindel2vcf, and skip_filter2
my %filter1_prefix;
my %filter2_prefix;
$filter1_prefix{'pass'} = "";
$filter1_prefix{'fail'} = "";
$filter2_prefix{'pass'} = "";
$filter2_prefix{'fail'} = "";
if ($paras{'apply_filter'} eq "true") { 
    $filter1_prefix{'pass'} = "CvgVafStrand_pass";
    $filter1_prefix{'fail'} = "CvgVafStrand_fail";
    $filter2_prefix{'pass'} = "Homopolymer_pass";
    $filter2_prefix{'fail'} = "Homopolymer_fail";
}

# run CvgVafStrandFilter
if ($paras{'apply_filter'} eq "true"  &&  $paras{'mode'} ne "pooled") {
    my $infn = "$var_file";
    my $passfn = "$var_file.$filter1_prefix{'pass'}";
    my $failfn = "$var_file.$filter1_prefix{'fail'}";
    if ( (not exists ($paras{'skip_filter1'})) or ($paras{'skip_filter1'} !~ /true/)) {
        print("Running filter 1 (CvgVafStrand).  Input: $infn  Pass Output: $passfn\n");
        CvgVafStrandFilter($infn, $passfn, $failfn, $paras{'mode'}, $paras{'min_coverages'}, 
            $paras{'min_var_allele_freq'}, ($paras{'require_balanced_reads'} =~ /true/), 
            ($paras{'remove_complex_indels'} =~ /true/), $paras{'child_var_allele_freq'}, 
            $paras{'parents_max_num_supporting_reads'}, $zero);
    } else {
        my $infn_base = basename($infn);
        print("Skipping filter 1 (CvgVafStrand).  Creating $passfn as link to $infn_base\n");
        my $return_code = system ( "ln -fs $infn_base $passfn" );
        die("Exiting ($return_code).\n") if $return_code != 0;
    }
}

# Conversion to VCF with pindel2vcf
my $skip_pindel2vcf = 0;
if ((not exists ($paras{'skip_pindel2vcf'})) or ($paras{'skip_pindel2vcf'} !~ /true/)) {
    my $var_src;
    my $var_dest;
    if ($paras{'apply_filter'} eq "true" && $paras{'mode'} ne "pooled") { # only case for filter1
        $var_src  = "$var_file.$filter1_prefix{'pass'}";
        $var_dest = "$var_file.$filter1_prefix{'pass'}.vcf";
    } else {
        $var_src  = "$var_file";
        $var_dest = "$var_file.vcf";
    }
    my $ref_base = `basename $paras{'REF'}`;
    chomp  $ref_base;
    my $logfile= ( exists($paras{'logfile'}) ? $paras{'logfile'} : "pindel2vcf.log" );

    print("Running pindel2vcf.pl  Input: $var_src  Output: $var_dest\n");
    runPindel2VCF($thisdir, $var_src, $var_dest, $ref_base, $logfile, $paras{'pindel2vcf'}, $paras{'REF'}, $paras{'date'},
        $paras{'heterozyg_min_var_allele_freq'}, $paras{'homozyg_min_var_allele_freq'});
} else {
    print("Skipping pindel2vcf\n");
}

# Filter 2: Homopolymer filtering
if ($paras{'apply_filter'} eq "true") {
    my $infn;
    my $passfn;
    my $failfn;

    if ($paras{'mode'} eq "pooled") {
        $infn = "$var_file.vcf";
        $passfn = "$var_file.$filter2_prefix{'pass'}.vcf";
        $failfn = "$var_file.$filter2_prefix{'fail'}.vcf";
    } else {
        $infn = "$var_file.$filter1_prefix{'pass'}.vcf";
        $passfn = "$var_file.$filter1_prefix{'pass'}.$filter2_prefix{'pass'}.vcf";
        $failfn = "$var_file.$filter1_prefix{'pass'}.$filter2_prefix{'fail'}.vcf";
    }

    if ((not exists ($paras{'skip_filter2'})) or ($paras{'skip_filter2'} !~ /true/)) {
        print("Running filter 2 (Homopolymer).  Input: $infn  Pass Output: $passfn\n");
        HomopolymerFilter($infn, $passfn, $failfn, $paras{'max_num_homopolymer_repeat_units'});
    } else {
        my $infn_base = basename($infn);
        print("Skipping filter 2 (Homopolymer).  Creating $passfn as link to $infn_base\n");
        my $return_code = system ( "ln -fs $infn_base $passfn" );
        die("Exiting ($return_code).\n") if $return_code != 0;
    }
}

1;
