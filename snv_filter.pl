#!/usr/bin/env perl
#----------------------------------
# @name GenomeVIP SNVs false-positives filter using VarScan
# @author Beifang Niu 
# @author R. Jay Mashl <rmashl@genome.wustl.edu>
# @author song cao
# @version 0.7 do not use tmp file: 
# @version 0.6 (rjm): make pass test explicit for extensibility
# @version 0.5 (rjm): workaround for stat on AWS
# @version 0.4 (rjm): handle zero-input case
# @version 0.3 (rjm): adjust filenames and parameter names; allow for commented lines
# @version 0.2 (rjm): added ENV and pass/fail splitter
# @version 0.1 (bn):  original
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
#use File::stat;


my $cmd;
my $input_fh_passed;
my $input_fh_failed;

# get paras from config file
my (%paras);
map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; print $v."\n";  $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<>);
# map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;
# print $paras{'output'}."\n";
#
print $paras{'variants_file'}."\n";
unless ( -e $paras{'variants_file'} ) { die "Variants input file could not be found !!! \n"; }

# Create readcounts input
my $input_fh = IO::File->new( $paras{'variants_file'} ) or die "Input variants file could not be opened.";
my $read_count_input  = $paras{'rc_in'};

my $read_count_input_fh = IO::File->new( $read_count_input, ">" ) or die "Temporary file could not be created. ";
my %seen = ();
map { unless( /^#/ ) { chomp; my @t = split /\t/; my $k = join( "\t", @t[0,1,1] ); $seen{$k} = 1;} } <$input_fh>;
map { $read_count_input_fh->print($_."\n") } keys %seen;
$input_fh->close        || die "Error on closing input variants file";

# Run readcounts
my $fs = `wc -l < $read_count_input`;
if( $fs !=  0)  {

	my $read_count_output =$paras{'rc_out'};

    #my ( undef, $read_count_output ) = tempfile();
    my $cmd_run_read_count = "$paras{'bam_readcount'} -w 10 -l $read_count_input -q $paras{'min_mapping_qual'} -b $paras{'min_base_qual'} -f $paras{'REF'}  $paras{'bam_file'}  > $read_count_output"; 
    print $cmd_run_read_count."\n";
    system( $cmd_run_read_count );
	my $fp_output_file = $paras{'fp_out'};

   # my ( undef, $fp_output_file ) = tempfile();
    my $cmd_run_varscan = "java $ENV{'JAVA_OPTS'} -jar $ENV{'VARSCAN_DIR'}/VarScan.jar fpfilter $paras{'variants_file'} $read_count_output --output-file $fp_output_file --keep-failures 1 --min-var-count $paras{'min_num_var_supporting_reads'} --min-var-freq $paras{'min_var_allele_freq'} --min-var-readpos $paras{'min_avg_rel_read_position'} --min-var-dist3 $paras{'min_avg_rel_dist_to_3prime_end'} --min-strandedness $paras{'min_var_strandedness'} --min-strand-reads $paras{'min_allele_depth_for_testing_strandedness'} --min-ref-basequal $paras{'min_ref_allele_avg_base_qual'} --min-var-basequal  $paras{'min_var_allele_avg_base_qual'} --max-rl-diff $paras{'max_rel_read_length_difference'} --max-var-mmqs $paras{'max_mismatch_qual_sum_for_var_reads'} --max-mmqs-diff $paras{'max_avg_mismatch_qual_sum_difference'} --min-ref-mapqual $paras{'min_ref_allele_avg_mapping_qual'} --min-var-mapqual $paras{'min_var_allele_avg_mapping_qual'} --max-mapqual-diff $paras{'max_avg_mapping_qual_difference'} ";
    print $cmd_run_varscan."\n";
    system( $cmd_run_varscan );
    
# (rjm) varscan --filtered-file option does not work at all, so manually extract failed calls
    $input_fh = IO::File->new( $fp_output_file  , "r")  or die "File for varscan fpfilter results could not be opened. ";
    $input_fh_passed = IO::File->new( $paras{'passfile'}, "w")  or die "could not open file for writing passed variants. ";
    $input_fh_failed = IO::File->new( $paras{'failfile'}, "w")  or die "could not open file for writing failed variants. ";
    while (<$input_fh>) {
	if( /^#/ ) {
	    $input_fh_passed->print($_);
	    $input_fh_failed->print($_);
	} else {
	    my @a = split /\t/;
	    if ( $a[6] eq "PASS" || $a[6] eq "." ) { # dot option is for gatk
		$input_fh_passed->print($_);
	    } else {
		$input_fh_failed->print($_);
	    }
	}
    }
    $input_fh_failed->close || die "Error on closing failed calls file";
    $input_fh_passed->close || die "Error on closing passed calls file";
    $input_fh->close        || die "Error on closing input variants file";
    #$cmd="rm -f $read_count_input $read_count_output $fp_output_file";
    #print $cmd."\n";
   # system( $cmd );

} else {

    print "# NOTICE: no variants are available for performing readcounts!\n";
    $input_fh = IO::File->new( $paras{'variants_file'} ) or die "Input variants file could not be opened.";
    $input_fh_passed = IO::File->new( $paras{'passfile'}, "w")  or die "could not open file for writing passed variants. ";
    $input_fh_failed = IO::File->new( $paras{'failfile'}, "w")  or die "could not open file for writing failed variants. ";
    while (<$input_fh>) {
	if( /^#/ ) {
	    $input_fh_passed->print($_);
	    $input_fh_failed->print($_);
	} 
    }
    $input_fh_failed->close || die "Error on closing failed calls file";
    $input_fh_passed->close || die "Error on closing passed calls file";
    $input_fh->close        || die "Error on closing input variants file";

}

    
1;
