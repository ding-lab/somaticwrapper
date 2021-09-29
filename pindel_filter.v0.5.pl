#!/usr/bin/env perl
#--------------------------------------
# @name GenomeVIP Pindel filters
# @author Beifang Niu
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# 
# @version 0.5 (rjm): improve documentation
# @version 0.4 (rjm): (optional) user-specified log file for appending; mode not required to be specified when not performing filtering
# @version 0.3 (rjm): generalize filter hierarchy filename handling and simply code structure
# @version 0.2 (rjm): added coverage, germline, trio, minimal pool filtering; revised approach; added failed pass. Adjusted filename and parameters names; allow for commented lines
# @version 0.1 (bn):  original somatic filter, written for (tumor,normal) column order
# updated by song cao to fix the bug about the somatic part#

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


#my $zero=0.001;
## song changes 0.001 to 0.01 ##
my $zero=0.01;

my $thisdir;
$thisdir=`dirname $ARGV[0]`;
chomp $thisdir;

# get paras from config file
my %paras; 
map { chomp; if ( !/^[#;]/ &&  /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
#map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;
if( $paras{'apply_filter'} eq "true" &&  !exists($paras{'mode'}) ) { die "Could not detect a filtering mode for filtering !!!\n"; }
# file exist ? 
unless ( -e $paras{'variants_file'} ) { die "input indels not exist !!! \n"; }


my $var_file        = $paras{'variants_file'};

# Filters for coverages, vaf, and balanced reads
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
my $filter1_pass_fh;
my $filter1_fail_fh;

# Conversion to VCF and homopolymer filtering
my $input_fh;
my $filter2_fh_pass;
my $filter2_fh_fail;



# Optional filter, part 1
if ($paras{'apply_filter'} eq "true"  &&  $paras{'mode'} ne "pooled") {
    $input_fh        = IO::File->new( "$var_file"                               ) or die "Could not open $var_file for reading $! ";
    $filter1_pass_fh = IO::File->new( "$var_file.$filter1_prefix{'pass'}",  ">" ) or die "Could not create $var_file.$filter1_prefix{'pass'} for writing $!";
    $filter1_fail_fh = IO::File->new( "$var_file.$filter1_prefix{'fail'}",  ">" ) or die "Could not create $var_file.$filter1_prefix{'fail'} for writing $!";


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
    if ($paras{'mode'} eq "germline") {
	while (<$input_fh>) {
	    chomp; 
	    my @t = split /\s+/;
	    if(  ($t[32] + $t[34] + $t[36] <  $paras{'min_coverages_n'}) || ($t[33] + $t[34] + $t[36] <  $paras{'min_coverages_n'})  ) {
		$filter1_fail_fh->print($_."\n");
		next;
	    }
	    if( ($t[34] + $t[36] + $t[34] + $t[36])/($t[32] + $t[33] + $t[34] + $t[36] + $t[34] + $t[36] ) <  $paras{'min_var_allele_freq'} ){
		$filter1_fail_fh->print($_."\n");
		next;
	    }
	    if($paras{'require_balanced_reads'} =~ /true/) {  
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
    # minimum coverage met for both samples; VAF threshold in tumor; zero variant support in normal; reads are considered balanced as long as there is nonzero read support in both directions in tumor
    if ($paras{'mode'} eq "somatic") {
	while (<$input_fh>) {
	    chomp; 
	    my @t = split /\s+/;
	    if(  ($t[32] + $t[34] + $t[36] <  $paras{'min_coverages_n'}) || ($t[33] + $t[34] + $t[36] <  $paras{'min_coverages_n'})  || ($t[39] + $t[41] + $t[43] < $paras{'min_coverages_t'}) ||  ($t[40] + $t[41] + $t[43] <  $paras{'min_coverages_t'})) {
		$filter1_fail_fh->print($_."\n");
		next;
	    }
### song changed ## 
       if( ($t[34] + $t[36] + $t[34] + $t[36])/($t[32] + $t[33] + $t[34] + $t[36] + $t[34] + $t[36] ) > $zero ||  ($t[41] + $t[43] + $t[41] + $t[43])/($t[39] + $t[40] + $t[41] + $t[43] + $t[41] + $t[43]) < $paras{'min_var_allele_freq'}) {

#	    if( ($t[34] + $t[36] + $t[34] + $t[36])/($t[32] + $t[33] + $t[34] + $t[36] + $t[34] + $t[36] ) <  $paras{'min_var_allele_freq'} ||  ($t[41] + $t[43] + $t[41] + $t[43])/($t[39] + $t[40] + $t[41] + $t[43] + $t[41] + $t[43]) > $zero) {
		$filter1_fail_fh->print($_."\n");
		next;
	    }
	    if($paras{'require_balanced_reads'} =~ /true/) {  
		if ( $t[41] == 0 || $t[43] == 0 ) {
		#if ( $t[34] == 0 ||  $t[36] == 0 || $t[41] > 0 || $t[43] > 0 ) {
		    $filter1_fail_fh->print($_."\n");
		    next;
		}
	    }
	    if($paras{'remove_complex_indels'} =~ /true/) {  
		if ( $t[1] eq "I" || $t[1] eq "D") {
		    if ( $t[1] eq "I" || ($t[1] eq "D" && $t[4] == 0) ) {
			print "Indel filter: passed\n";
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
    # minimum coverage met for all samples; VAF threshold in child; maximum allowed variant support in parents combined; reads are considered balanced as long as there is nonzero read support in both directions in child

  ## not needed for somatic variant callings ##

   if ($paras{'mode'} eq "trio") {   # trio
	while (<$input_fh>) {
	    chomp; 
	    my @t = split /\s+/;
	    if(  ($t[32] + $t[34] + $t[36] <  $paras{'min_coverages_n'}) || ($t[33] + $t[34] + $t[36] <  $paras{'min_coverages_n'})  || ($t[39] + $t[41] + $t[43] < $paras{'min_coverages_n'}) ||  ($t[40] + $t[41] + $t[43] <  $paras{'min_coverages_n'}) ||   ($t[46] + $t[48] + $t[50] < $paras{'min_coverages_n'}) ||  ($t[47] + $t[48] + $t[50] <  $paras{'min_coverages_n'}) ) {
		$filter1_fail_fh->print($_."\n");
		next;
	    }
	    if( ($t[48] + $t[50] + $t[48] + $t[50])/($t[46] + $t[47] + $t[48] + $t[50] + $t[48] + $t[50]) < $paras{'child_var_allele_freq'}    ) { 
		$filter1_fail_fh->print($_."\n");
		next;
	    }
	    if( ($t[34] + $t[36]) + ($t[41] + $t[43]) >  $paras{'parents_max_num_supporting_reads'}  ) {
		$filter1_fail_fh->print($_."\n");
		next;
	    }
	    if($paras{'require_balanced_reads'} =~ /true/) {  
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




# run pindel2vcf
my $pindel2vcf_command = "";
my $ref_base = `basename $paras{'REF'}`;
chomp  $ref_base;
my $result;
my $var_src;
my $var_dest;

if ($paras{'apply_filter'} eq "true" && $paras{'mode'} ne "pooled") { # only case for filter1
    $var_src  = "$var_file.$filter1_prefix{'pass'}";
    $var_dest = "$var_file.$filter1_prefix{'pass'}.vcf";
} else {
    $var_src  = "$var_file";
    $var_dest = "$var_file.vcf";
}
my $logfile="pindel2vcf.log";
if(exists($paras{'logfile'})) { $logfile=$paras{'logfile'}; }
# NOTE: pindel -co option seems to work on vcf output but warning messages may happend regardless
$pindel2vcf_command = "$paras{'pindel2vcf'} -r $paras{'REF'} -R $ref_base -p $var_src -d $paras{'date'} -v $var_dest -he $paras{'heterozyg_min_var_allele_freq'} -ho $paras{'homozyg_min_var_allele_freq'}  >> $thisdir/$logfile 2>&1";
# print $pindel2vcf_command."\n";
$result = system( $pindel2vcf_command );

    
# Optional filter, part 2; here pooled is ok
if ($paras{'apply_filter'} eq "true") {
    
    if ($paras{'mode'} eq "pooled") {
	$input_fh        = IO::File->new( "$var_file.vcf"                             ) or die "Could not open $var_file.vcf for reading $!";
	$filter2_fh_pass = IO::File->new( "$var_file.$filter2_prefix{'pass'}.vcf", ">") or die "Could not create $var_file.$filter2_prefix{'pass'}.vcf for writing $!";
	$filter2_fh_fail = IO::File->new( "$var_file.$filter2_prefix{'fail'}.vcf", ">") or die "Could not create $var_file.$filter2_prefix{'fail'}.vcf for writing $!";
    } else {
	$input_fh        = IO::File->new( "$var_file.$filter1_prefix{'pass'}.vcf"                             ) or die "Could not open $var_file.$filter1_prefix{'pass'}.vcf for reading $!";
	$filter2_fh_pass = IO::File->new( "$var_file.$filter1_prefix{'pass'}.$filter2_prefix{'pass'}.vcf", ">") or die "Could not create $var_file.$filter1_prefix{'pass'}.$filter2_prefix{'pass'}.vcf for writing $!";
	$filter2_fh_fail = IO::File->new( "$var_file.$filter1_prefix{'pass'}.$filter2_prefix{'fail'}.vcf", ">") or die "Could not create $var_file.$filter1_prefix{'pass'}.$filter2_prefix{'fail'}.vcf for writing $!";
    }
    
    while ( <$input_fh> ) {
#	print;
	if ( /^#/ ) { $filter2_fh_pass->print($_); next };
	my @a= split /\t/; 
	my @b = split/\;/, $a[7]; 
	for ( my $i=0; $i<scalar(@b); $i++) { 
	    if ( $b[$i]=~/^HOMLEN/ ) { 
		my @c = split/=/, $b[$i]; 
		if ( $c[1] <= $paras{'max_num_homopolymer_repeat_units'} ) { 
		    $filter2_fh_pass->print($_); 
		}  else {
		    $filter2_fh_fail->print($_); 
		}
		last;
	    } 
	}
    }
    $filter2_fh_fail->close;
    $filter2_fh_pass->close;
    $input_fh->close;
    
}

1;
