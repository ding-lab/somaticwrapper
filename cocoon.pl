#!/usr/bin/perl

use strict;
use warnings;


use Getopt::Long;
my $mapq = 10;			#by default
my $distance = 5;		#by default
my $freq_cooccur = 0.8; 	#by default
my $vaf_diff = 10;		#by default
my $mergeFlag = 0;		#by default
my $snvonlyFlag = 0;		#by default
my $helpFlag = 0;
my $genome;
my $gtf;
my $bam;
my $samtools;

GetOptions(	"distance=n" 		=> \$distance,
		"mapq=n"		=> \$mapq,
		"samt=s"  => \$samtools,
		"freq_cooccur=n" 	=> \$freq_cooccur,
		"vaf_diff=n"		=> \$vaf_diff,
		"merge"			=> \$mergeFlag,
		"bam=s"			=> \$bam,
		"genome=s"		=> \$genome,
		"gtf=s"			=> \$gtf,
		"snvonly"		=> \$snvonlyFlag,
		"help" 			=> \$helpFlag );

if(!defined($ARGV[2]) || $helpFlag)
{
	print "\nUsage: perl cocoons.pl MAF OUTPUT DIRLOG <OPTIONS>

COrrecting CO-Occuring multi-Nucleotide variationS (COCOONS) in Mutation Annotation Format (MAF) file

REQUIRED:
	MAF			MAF file to be inspected
	OUTPUT			Output file name (used as prefix for merged file in MERGING mode)
	DIRLOG 		 Directory for outputing GTF_CDS and CDS_BED files 
OPTIONS:
	--distance		Maximum distance between two mutations to be considered [default: 5]
	--mapq			Minimum read mapping quality; used if given bam [default: 10]
	--snvonly		Only SNVs are considered [default: no]
	--freq_cooccur		Minimum cooccuring frequency of the adjacent mutations; used if given bam [default: 0.8]
	--vaf_diff		Maximum variant allele frequency (VAF) difference to define DNP; used if not given bam [default:10]
	--bam			File containing bam file locations for the samples used in MAF (format:SAMPLE<TAB>BAM_LOCATION in each line)
	--merge			Merge the adjacent mutations (MERGING mode) or not [default: no]
	--genome		Human genome sequence (fasta) [required in MERGING mode]
	--gtf			Ensembl human annotation (gtf) [required in MERGING mode]
	--help			Print this help message

*** Questions \& Bug Reports: Qingsong Gao (qingsong.gao\@wustl.edu)

";
	exit
}


##########################################
########### checking input files #########
#############=> START <=##################

my $maf = $ARGV[0];
my $output =$ARGV[1];
my $dirlog=$ARGV[2];

 
my $f_cds_bed; 
my $f_gtf_cds; 
open(OUT, ">$output");

open(MAF, "$maf");
my $header = <MAF>;
print OUT "$header";

my $t_depth_idx = -1;
my $t_alt_count_idx = -1;

chomp $header;
my @h = split(/\t/, $header);
for(my $i=0;$i<=$#h; $i++)
{
        if($h[$i] eq "t_depth")
        {
                $t_depth_idx = $i;
        }elsif($h[$i] eq "t_alt_count")
        {
                $t_alt_count_idx = $i;
        }
}

my %bam_location;

if(!$bam)
{
#	print STDERR "BAM files not provided\n";
        if($t_depth_idx<0 || $t_alt_count_idx<0)
        {
                die "If bam files are not available, \"t_depth\" and \"t_alt_count\" columns in MAF are required!\n";
        }
}else
{
#	print STDERR "BAM file provided\n";
	open(BAM, "$bam");
	while(<BAM>)
	{
		chomp;
		my @line = split(/\s+/, );
		$bam_location{$line[0]} = $line[1];
	}
}

my %cdna;
my %direct;
if($mergeFlag)
{
        if($genome && $gtf)
        {
                open(OUTPUT, ">$output.merge");
		print OUTPUT "$header";

		### GTF to CDS GTF ###
		my %left;
		my %right;
		my ($tmp, $tmp1, $tmp2);
		open(GTF, "$gtf");
		 $f_gtf_cds=$dirlog."/GTF_CDS";
		open(CDSOUT, ">$f_gtf_cds");

		while(<GTF>)
		{
	        	chomp;
        		next if $_ =~m/^#/;
        		my @l=split(/\t/,);
        		my @a=split(/\"/,$l[8]);
         ### $a[3], transcript id; $a[5], exon number is a[9], not a[5] ##
                        #print $a[5],"\t",$a[9],"\n";
                        #<STDIN>;

	        	if($l[2] eq "CDS" || $l[2] eq "stop_codon")
        		{
                        #print $a[5],"\t",$a[9],"\n";
                        #<STDIN>;
                		if(! exists $left{$a[5]}{$a[9]})
                		{
                		        $left{$a[5]}{$a[9]}=$l[3];
                		}else
                		{
                		        if($l[3]<$left{$a[5]}{$a[9]})
                		        {
                		                $left{$a[5]}{$a[9]}=$l[3];
                		        }
                		}
                		if(! exists $right{$a[5]}{$a[9]})
                		{
                		        $right{$a[5]}{$a[9]}=$l[4];
                		}else
                		{
                		        if($l[4]>$right{$a[5]}{$a[9]})
                		        {
                		                $right{$a[5]}{$a[9]}=$l[4];
                		        }
                		}
        		}
		}
		open(GTF, "$gtf");
                while(<GTF>)
                {
                        chomp;
                        next if $_ =~m/^#/;
                        my @l=split(/\t/,);
                        my @a=split(/\"/,$l[8]);
                        next if $l[2] ne "exon";
                        next if !exists $left{$a[5]};
                        next if !exists $left{$a[5]}{$a[9]};

                        #print $a[5],"\t",$a[9],"\n";
                        #<STDIN>;

                        if($l[6] eq "+")
                        {
                                if($l[3]==$left{$a[5]}{$a[9]})
                                {
                                        if($right{$a[5]}{$a[9]}==$l[4])
                                        {
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$l[3]\t$l[4]\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }else
                                        {
                                                $tmp=$right{$a[5]}{$a[9]}+1;
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$l[3]\t$right{$a[5]}{$a[9]}\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }
                                }else
                                {
                                        if($right{$a[5]}{$a[9]}==$l[4])
                                        {
                                                $tmp=$left{$a[5]}{$a[9]}-1;
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$left{$a[5]}{$a[9]}\t$l[4]\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }else
                                        {
                                                $tmp1=$left{$a[5]}{$a[9]}-1;
                                                $tmp2=$right{$a[5]}{$a[9]}+1;
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$left{$a[5]}{$a[9]}\t$right{$a[5]}{$a[9]}\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }
                                }
                        }else    
                        {
                                if($l[3]==$left{$a[5]}{$a[9]})
                                {
                                        if($right{$a[5]}{$a[9]}==$l[4])
                                        {
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$l[3]\t$l[4]\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }else
                                        {
                                                $tmp=$right{$a[5]}{$a[9]}+1;
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$l[3]\t$right{$a[5]}{$a[9]}\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }
                                 }else
                                 {
                                        if($right{$a[5]}{$a[9]}==$l[4])
                                        {
                                                $tmp=$left{$a[5]}{$a[9]}-1;
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$left{$a[5]}{$a[9]}\t$l[4]\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }else
                                        {
                                                $tmp1=$left{$a[5]}{$a[9]}-1;
                                                $tmp2=$right{$a[5]}{$a[9]}+1;
                                                print CDSOUT "$l[0]\t$l[1]\tCDS\t$left{$a[5]}{$a[9]}\t$right{$a[5]}{$a[9]}\t$l[5]\t$l[6]\t$l[7]\t$l[8]\n";
                                        }
                                }
                        }
                }
                close CDSOUT;
		### GTF to BED ###
		open(GTFCDS, "$f_gtf_cds");

		$f_cds_bed=$dirlog."/CDS_BED"; 

		open(OUTBED, ">$f_cds_bed");

		my %block_count;
		my %strand;
		my %starts;
		my %ends;
		my %chr;
		while(<GTFCDS>)
		{
			chomp;
			my @line = split(/\t/,);
			my @ids = split /\"/, $line[8];
					## hg38 refernece ##
                	my $trans_id = $ids[5];
					#print $trans_id,"\n";
					#<STDIN>;
                	$block_count{$trans_id}++;
                	my $tmp_start = $line[3]-1;
                	my $tmp_end = $line[4];
                	$strand{$trans_id}=$line[6];
                	$starts{$trans_id}.="$tmp_start ";
                	$ends{$trans_id}.="$tmp_end ";
                	$chr{$trans_id} = "$line[0]";
		}
		foreach my $id (keys %block_count)
		{
        		my @lens = ();
        		my @block_start = ();
        		my @start = split(/ /, $starts{$id});
        		my @end = split(/ /, $ends{$id});
        		my @sorted_start = sort {$a <=> $b} @start;
        		my @sorted_end = sort {$a <=> $b} @end;
       		 	for(my $i=0; $i<=$#start; $i++)
        		{
                		my $len = $sorted_end[$i]-$sorted_start[$i];
                		push @lens, $len;
                		my $re_start = $sorted_start[$i]-$sorted_start[0];
                		push @block_start, $re_start;
        		}
        		my $size = join ",", @lens;
        		my $block_starts = join ",", @block_start;
        		print OUTBED "$chr{$id}\t$sorted_start[0]\t$sorted_end[-1]\t$id\t0\t$strand{$id}\t$sorted_start[0]\t$sorted_start[0]\t0\t$block_count{$id}\t$size\t$block_starts\n";
		}
		close OUTBED;

		#`rm $f_cds_bed`; 
		#`rm $f_gtf_cds`;

		### BED to FA ###

		open(REF,"$genome");
		my %human;
		my $id;

		while(<REF>)
		{
        		chomp;
        		my $first = substr($_,0,1);
        		if($first eq ">")
        		{
                		my @line = split(/\s+/, substr($_,1));
                		$id=$line[0];
						$id=~s/chr//g;  ### remove chr ##
        		}else
        		{
                		$human{$id}.=uc($_);
        		}
		}
		close REF;
		open(DA,"<$f_cds_bed");
		while(<DA>)
		{
        		chomp;
        		my @l=split(/\t/,);
        		my @size=split(/\,/,$l[10]);
        		my @start=split(/\,/,$l[11]);
        		my $seq="";

            	#print $l[3],"\t",$seq,"\n";
				#<STDIN>;
        		next if !exists $human{$l[0]};
        		for(my $i=0;$i<=$#size;$i++)
        		{
                		$seq.=substr($human{$l[0]},$l[1]+$start[$i],$size[$i]);
        		}
        		$seq=uc($seq);
        		if($l[5] eq "-")
        		{
                		$seq=~tr/[AGCT]/[TCGA]/;
                		$seq=reverse($seq);
        		}
			$direct{$l[3]} = $l[5];
			$cdna{$l[3]} = $seq;
			#print $l[3],"\t",$seq,"\n";
			#<STDIN>;
		}
        }else
        {
                die "Genome sequence(fasta) and annotation(gtf) are required in MERGE mode!\n";
        }
		close DA;
		 #`rm $f_cds_bed`;
        #`rm $f_gtf_cds`;
}

###############=> END <=##################
##########################################
#foreach $a (keys %human) {print "$a\n";}

##########################################
## fuctions to decipher variants in BAM ##
#############=> START <=##################

sub to_array
{
        my $input=shift;
        my @out=();
        my $num=0;
        my @arr=split(//,$input);
        for(my $i=0;$i<=$#arr;$i++)
        {
                if($arr[$i]=~m/[0-9]/)
                {
                        $num=$num*10+$arr[$i];  
                }else
                {
                        push @out,$num;
                        push @out,$arr[$i];
                        $num=0;
                }
        }
        my $last=substr($input,-1);
        if($last=~m/[0-9]/)
        {
                push @out,$num;
        }
        return @out;
}

sub to_arrays
{
        my $input=shift;
        if($input!~m/\^/)
        {
                return to_array($input);
        }else
        #with deletion
        {
                my @out=();
                my @array=($input=~m/(\^\D+)/g);
                for(my $i=0;$i<=$#array;$i++)
                {
                        my @sep_array=split(/\Q$array[$i]/,$input);
                        my $num=scalar @sep_array;
                        $input=join("$array[$i]",@sep_array[1..$#sep_array]);   
			#to avoid same deletions
                        @out=(@out,to_array($sep_array[0]));
                        my $del=substr($array[$i],1);
                        push @out,"\^".$del;
                }
                @out=(@out,to_array($input));
        }
}

sub misPos
{
        my @array=@_;
        for(my $i=10;$i<=$#array;$i++)
        {
                if($array[$i]=~m/MD:Z:/)
                {
                        my @md=split(/\:/,$array[$i]);
                        return $md[2];
                        last;   
                }
        }
}

sub reports
{
	my $align = shift;
	my @line = split(/\t/, $align);
	my $seq = $line[9];
	my $start = $line[3];
	my %variants;
	my @out;

	my $cigar = $line[5];
        my @number = ($cigar=~m/(\d+)\w/g);
        my @type = ($cigar=~m/\d+(\w)/g);
	my $tmp_cigar = "";

	if($cigar=~m/S/)
	#soft-clipped
	{
		if($type[-1] eq "S")
		#last section goes first if softclipped happen on both ends;
		{
			$seq=substr($seq, 0, length($seq)-$number[-1]);
			for(my $k=0; $k<$#type; $k++)
			{
				$tmp_cigar.=$number[$k].$type[$k];
			}
			pop @type;
			pop @number;
			$cigar = $tmp_cigar;
		}
		$tmp_cigar ="";
		if($type[0] eq "S")
		{
			$seq=substr($seq, $number[0]);
			for(my $a=1; $a<=$#type; $a++)
			{
				$tmp_cigar.=$number[$a].$type[$a];
			}
			shift @type;
			shift @number;
			$cigar = $tmp_cigar;
		}
		die "Softclip can only happen at very end!" if $cigar=~m/S/;
	}

	#decipher cigar string to report insertions
	my $cigar_pos = 0;
        my $adj_de = 0;
        if($cigar=~m/I/)
        {
                my $new_seq = "";
                for(my $i=0; $i<=$#type; $i++)
                {
                        if($type[$i] eq "I")
                        {
				$variants{$start+$cigar_pos+$adj_de-1} = "\-\>".substr($seq, $cigar_pos, $number[$i]);
				$cigar_pos += $number[$i];
                        }elsif($type[$i] eq "D")
                        {
                                $adj_de += $number[$i];
                        }else
                        {
                                $new_seq.=substr($seq, $cigar_pos, $number[$i]);
                                $cigar_pos += $number[$i];
                        }
                }
                $seq = $new_seq;
        }
	#decipher MD string to report deletion and mismatch
	my $md = misPos(@line);
	my @md_array =  to_arrays($md);
	my $md_pos = 0;
	my $adj_del = 0;
	for(my $j=0; $j<=$#md_array; $j++)
	{
		if($md_array[$j]=~m/\d/)
		{
			$md_pos += $md_array[$j];
		}else
		{
			if($md_array[$j]=~m/\^/)
			{
				$variants{$start+$md_pos} = substr($md_array[$j],1)."\>\-";
				$adj_del += length(substr($md_array[$j],1));
				$md_pos += length(substr($md_array[$j],1));
			}else
			{
				$variants{$start+$md_pos} = $md_array[$j]."\>".substr($seq, $md_pos-$adj_del, 1);
				$md_pos++;
			}
		}
	}

	foreach my $p (sort {$a<=>$b} keys %variants)
	{
		push @out,"$p\:$variants{$p}\;";
	}
	return join("\t", @out);
}

###############=> END <=##################
##########################################

##########################################
## fuctions to compare reads in MAF/BAM ##
##############=> START <=##################

sub readsConfirmed
{
	my $input = shift;
	my @line = split(/\n/, $input);
	my ($first, $middle, $last, @f, @m, @l, $seq, @sequences, @reads, $mismatch, $var1, $var2, $var3);
	my $and = 0;
	my $or = 0.000001;
	my $realD;
	my $realN;

	if(scalar @line > 3)
	{
#		print STDERR "Skipped because >3 variants cannot affect the same codon\n";
#		print STDERR "To include the ones on the same codon from this cluster, set \"--distance 2\"\n";
#		print STDERR "This could be taken care of in next version if necessary\n";
#		print STDERR "The minimum distance between any two variation in this cluster will be used\n";
	}elsif(scalar @line == 3)
	{
                $first = $line[0];
                $middle = $line[1];
                $last = $line[2];
                @f = split(/\t/, $first);
                @m = split(/\t/, $middle);
                @l = split(/\t/, $last);
                if(!$bam)
                {
                        if(abs($f[$t_alt_count_idx]/$f[$t_depth_idx]-$m[$t_alt_count_idx]/$m[$t_depth_idx]) <= $vaf_diff && abs($l[$t_alt_count_idx]/$l[$t_depth_idx]-$m[$t_alt_count_idx]/$m[$t_depth_idx]) <= $vaf_diff)
                        {
				$realN = 3;
				$realD = $l[6] - $f[5];
                                return "\#\t$realN\t$realD\n".$input."\n";
                        }else
                        {
                                if(abs($f[$t_alt_count_idx]/$f[$t_depth_idx]-$m[$t_alt_count_idx]/$m[$t_depth_idx]) <= $vaf_diff)
                                {
					$realN = 2;
					$realD = $m[6] - $f[5];
                                        return "\#\t$realN\t$realD\n".join("\n", @line[0..1])."\n";
                                }elsif(abs($l[$t_alt_count_idx]/$l[$t_depth_idx]-$m[$t_alt_count_idx]/$m[$t_depth_idx]) <= $vaf_diff)
                                {
					$realN = 2;
					$realD = $l[6] - $m[5];
                                        return "\#\t$realN\t$realD\n".join("\n", @line[1..2])."\n";
                                }elsif(abs($f[$t_alt_count_idx]/$f[$t_depth_idx]-$l[$t_alt_count_idx]/$l[$t_depth_idx]) <= $vaf_diff)
				{
					$realN = 2;
					$realD = $l[6] - $f[5];
					return "\#\t$realN\t$realD\n".join("\n", $line[0], $line[2])."\n";
				}
                        }
                }else
		{
                        $seq = `$samtools view -q $mapq $bam_location{$f[15]} $f[4]:$f[5]-$l[6]`;
                        @sequences = split(/\n/, $seq);
                        $var1 = $f[5].":".$f[10]."\>".$f[12];
                        $var2 = $m[5].":".$m[10]."\>".$m[12];
                        $var3 = $l[5].":".$l[10]."\>".$l[12];
                        my $and12 = 0;
                        my $and23 = 0;
			my $and13 = 0;
                        my $or12 = 0.000001;
                        my $or23 = 0.000001;
			my $or13 = 0.000001;
                        foreach my $s (@sequences)
                        {
                                @reads = split(/\t/,$s);
				next if $reads[5]=~m/H/;
                                if($reads[3]<=$f[5] && $reads[3]+length($reads[9])-1>=$f[5] && $reads[3]<=$l[5] && $reads[3]+length($reads[9])-1>=$l[6])
				{
                                        $mismatch = reports($s);
                                        if($mismatch=~m/$var1/ && $mismatch=~m/$var2/ && $mismatch=~m/$var3/)
                                        {
                                                $and++;
                                        }
                                        if($mismatch=~m/$var1/ || $mismatch=~m/$var2/ || $mismatch=~m/$var3/)
                                        {
                                                $or++;
                                        }
                                        if($mismatch=~m/$var1/ && $mismatch=~m/$var2/)
                                        {
                                                $and12++;
                                        }
                                        if($mismatch=~m/$var2/ && $mismatch=~m/$var3/)
                                        {
                                                $and23++;
                                        }
                                        if($mismatch=~m/$var1/ && $mismatch=~m/$var3/)
                                        {
                                                $and13++;
                                        }
                                        if($mismatch=~m/$var1/ || $mismatch=~m/$var2/)
                                        {
                                                $or12++;
                                        }
                                        if($mismatch=~m/$var2/ || $mismatch=~m/$var3/)
                                        {
                                                $or23++;
                                        }
                                        if($mismatch=~m/$var1/ || $mismatch=~m/$var3/)
                                        {
                                                $or13++;
                                        }
				}
			}
                        if($and/$or>=$freq_cooccur && $and>=2)
                        {
                                $realN = 3;
                                $realD = $l[6] - $f[5];
                                return "\#\t$realN\t$realD\n".$input."\n";
                        }else
                        {
                                if($and12/$or12>=$freq_cooccur && $and12>=2)
                                {
                                        $realN = 2;
                                        $realD = $m[6] - $f[5];
                                        return "\#\t$realN\t$realD\n".join("\n", @line[0..1])."\n";
                                }elsif($and23/$or23>=$freq_cooccur && $and23>=2)
                                {
                                        $realN = 2;
                                        $realD = $l[6] - $m[5];
                                        return "\#\t$realN\t$realD\n".join("\n", @line[1..2])."\n";
                                }elsif($and13/$or13>=$freq_cooccur && $and13>=2)
				{
                                        $realN = 2;
                                        $realD = $l[6] - $f[5];
					return "\#\t$realN\t$realD\n".join("\n", $line[0], $line[2])."\n";
				}
                        }
		}
	}elsif(scalar @line == 2)
	{
		$first = $line[0];
		$last = $line[1];
		@f = split(/\t/, $first);
		@l = split(/\t/, $last);

		if(!$bam)
		{
			if(abs($f[$t_alt_count_idx]/$f[$t_depth_idx]-$l[$t_alt_count_idx]/$l[$t_depth_idx]) <= $vaf_diff)
			{
				$realN = 2;
				$realD = $l[6] - $f[5];
				return "\#\t$realN\t$realD\n".$input."\n";
			}
		}else
		{
			$seq = `$samtools view -q $mapq $bam_location{$f[15]} $f[4]:$f[5]-$l[6]`;
			@sequences = split(/\n/, $seq);
			$var1 = $f[5].":".$f[10]."\>".$f[12];
			$var2 = $l[5].":".$l[10]."\>".$l[12];
			#print $var1,"\t",$var2,"\n"; 
#			<STDIN>;	
			foreach my $t (@sequences)
			{
				#print $t,"\n"; 
				#<STDIN>; 
				@reads = split(/\t/,$t);
				next if $reads[5]=~m/H/;
        			#ignore hard clipped reads
				if($reads[3]<=$f[5] && $reads[3]+length($reads[9])-1>=$f[5] && $reads[3]<=$l[5] && $reads[3]+length($reads[9])-1>=$l[6])
				#samtools get the reads overlapping with given region
				#this is to ensure the reads cover both sites
				{
					$mismatch = reports($t);
					if($mismatch=~m/$var1/ && $mismatch=~m/$var2/)
					{
						$and++;
					}	
					if($mismatch=~m/$var1/ || $mismatch=~m/$var2/)
					{
						$or++;
					}
				}
			}
			#print $and,"\t",$or,"\n"; 
			#<STDIN>;
			if($and/$or>=$freq_cooccur && $and>=2)
			{
				$realN = 2;
                                $realD = $l[6] - $f[5];
                                return "\#\t$realN\t$realD\n".$input."\n";
			}
		}
	}else
	{
		die "There should be at least two mutations\n";
	}
}

###############=> END <=##################
###########################################

##########################################
## fuctions for merging and annotation ###
#############=> START <=##################

my %genetic_code = (
        'TTT'=>'F', #Phenylalanine
        'TTC'=>'F', #Phenylalanine
        'TTA'=>'L', #Leucine
        'TTG'=>'L', #Leucine
        'TCA'=>'S', #Serine
        'TCC'=>'S', #Serine
        'TCG'=>'S', #Serine
        'TCT'=>'S', #Serine
        'TAC'=>'Y', #Tyrosine
        'TAT'=>'Y', #Tyrosine
        'TAA'=>'*', #Stop
        'TAG'=>'*', #Stop
        'TGC'=>'C', #Cysteine
        'TGT'=>'C', #Cysteine
        'TGA'=>'*', #Stop
        'TGG'=>'W', #Tryptophan

        'CTA'=>'L', #Leucine
        'CTC'=>'L', #Leucine
        'CTG'=>'L', #Leucine
        'CTT'=>'L', #Leucine
        'CCA'=>'P', #Proline
        'CCC'=>'P', #Proline
        'CCG'=>'P', #Proline
        'CCT'=>'P', #Proline
        'CAT'=>'H', #Histidine
        'CAC'=>'H', #Histidine
        'CAA'=>'Q', #Glutamine
        'CAG'=>'Q', #Glutamine
        'CGA'=>'R', #Arginine
        'CGC'=>'R', #Arginine
        'CGG'=>'R', #Arginine
        'CGT'=>'R', #Arginine

        'ATA'=>'I', #Isoleucine
        'ATC'=>'I', #Isoleucine
        'ATT'=>'I', #Isoleucine
        'ATG'=>'M', #Methionine
        'ACA'=>'T', #Threonine
        'ACC'=>'T', #Threonine
        'ACG'=>'T', #Threonine
        'ACT'=>'T', #Threonine
        'AAC'=>'N', #Asparagine
        'AAT'=>'N', #Asparagine
        'AAA'=>'K', #Lysine
        'AAG'=>'K', #Lysine
        'AGC'=>'S', #Serine#Valine
        'AGT'=>'S', #Serine
        'AGA'=>'R', #Arginine
        'AGG'=>'R', #Arginine

        'GTA'=>'V', #Valine
        'GTC'=>'V', #Valine
        'GTG'=>'V', #Valine
        'GTT'=>'V', #Valine
        'GCA'=>'A', #Alanine
        'GCC'=>'A', #Alanine
        'GCG'=>'A', #Alanine
        'GCT'=>'A', #Alanine
        'GAC'=>'D', #Aspartic Acid
        'GAT'=>'D', #Aspartic Acid
        'GAA'=>'E', #Glutamic Acid
        'GAG'=>'E', #Glutamic Acid
        'GGA'=>'G', #Glycine
        'GGC'=>'G', #Glycine
        'GGG'=>'G', #Glycine
        'GGT'=>'G', #Glycine
 );

my %codon_abbr = (
	'F'=>'Phe',
	'L'=>'Leu',
	'S'=>'Ser',
	'Y'=>'Tyr',
	'X'=>'Ter',
	'C'=>'Cys',
	'W'=>'Trp',
	'P'=>'Pro',
	'H'=>'His',
	'Q'=>'Gln',
	'R'=>'Arg',
	'I'=>'Ile',
	'M'=>'Met',
	'T'=>'Thr',
	'N'=>'Asn',
	'K'=>'Lys',
	'V'=>'Val',
	'A'=>'Ala',
	'D'=>'Asp',
	'E'=>'Glu',
	'G'=>'Gly',
        '*'=>'Ter',
);

sub uniq
{
        my %seen;
	return grep { !$seen{$_}++ } @_;
}

sub toMerge
{
        my $input = shift;
        my @line = split(/\n/, $input);
	my @h = split(/\t/,$line[0]);
	my $result;
	my @trans;
	my @amin;
	my ($p, $a, $poa, $amino, $ref_dna, $ref_aa, $ref_given, $effect, $alt_aa, $alt_dna);
	my $new_start = 1000000000;
	my $new_end = 0;
	my $seq = "";
	my @l;
	if($h[2] > 2)
	{
		$result=join("\t", @h)."\tMerged\=No\(far away\)\n".join("\n", @line[1..$#line])."\n";
		return $result;
	}else
	{
		shift @line;
		if($input=~m/Splice_Site/ || $input=~m/Translation_Start_Site/ || $input=~m/Nonstop_Mutation/)
		{
			$result=join("\t", @h)."\tMerged\=No\(not cds\)\n".join("\n", @line[0..$#line])."\n";
			return $result;
		}else
		{
			foreach my $record (@line)
			{
				my @l = split(/\t/,$record);
				push @trans, $l[37];
				if($l[36]=~m/p.\D+(\d+)\D+/)
				{
					push @amin,$1;
				}else
				{
					push @amin,rand(1);
				}
			}
			my @uniq_trans = uniq(@trans);
			my @uniq_amin = uniq(@amin);
			if(scalar @uniq_trans >1)
			{
				$result=join("\t", @h)."\tMerged\=No\(different transcripts\)\n".join("\n", @line[0..$#line])."\n";
				return $result;
			}else
			{
				if(scalar @uniq_amin > 1)
				{
					$result=join("\t", @h)."\tMerged\=No\(different codons\)\n".join("\n", @line[0..$#line])."\n";
					return $result;
				}else
				{
					foreach my $record (@line)
					{
						@l = split(/\t/,$record);
						$new_start = ($new_start > $l[5]) ? $l[5] : $new_start;
						$new_end = ($new_end < $l[6]) ? $l[6] : $new_end;
						$p = substr($l[34], 2, -3);
						$a = substr($l[34], -1);
						$poa = substr($l[36], 3, -1);
						$amino = substr($l[36], 3, -1)*3-3;
						if (substr($l[36], 3, -1)*3<$p || substr($l[36], 3, -1)*3-2>$p)
						{
							$result=join("\t", @h)."\tMerged\=No\(wrong annotation in maf\)\n".join("\n", @line[0..$#line])."\n";
							return $result;
						}else
						{
							if(!exists $cdna{$l[37]})
							{
								$result=join("\t", @h)."\tMerged\=No\(wrong annotation in maf, not exist cdna, $l[37]\)\n".join("\n", @line[0..$#line])."\n";
								return $result;
							}
							$ref_dna = substr($cdna{$l[37]}, $amino, 3);
							if(length($ref_dna) !=3)
							{
								$result=join("\t", @h)."\tMerged\=No\(wrong annotation in maf, ref_dna not equal to 3 \)\n".join("\n", @line[0..$#line])."\n";
								return $result;
							}else
							{
								$ref_aa=$genetic_code{$ref_dna};
								$ref_given=substr($l[36], 2, 1);
								if($ref_aa ne $ref_given)
								{
									$result=join("\t", @h)."\tMerged\=No\(wrong annotation in maf, ref_aa not equal to ref_given \)\n".join("\n", @line[0..$#line])."\n";
									return $result;
								}else
								{
									if($seq eq "")
									{
										$seq = $cdna{$l[37]};
										substr($seq, $p-1, 1)=$a;
									}else
									{
										substr($seq, $p-1, 1)=$a;
									}
								}
							}
						}
					}
					$alt_dna = substr($seq, $amino, 3);
					$alt_aa = $genetic_code{$alt_dna};
					if($ref_aa eq $alt_aa)
					{
						$effect = "Silent";
					}else
					{
						if($alt_aa eq "\*")
						{
							$effect = "Nonsense_Mutation";
						}else
						{
							$effect = "Missense_Mutation";
						}
					}
					if($direct{$l[37]} eq "-")
					{
						$ref_dna=~tr/[AGCT]/[TCGA]/;
						$ref_dna=reverse($ref_dna);
						$alt_dna=~tr/[AGCT]/[TCGA]/;
						$alt_dna=reverse($alt_dna);
					}
                                        if($h[2]==1)
                                        {
                                                if(substr($ref_dna, 0, 1) eq substr($alt_dna, 0, 1))
                                                {
                                                        $ref_dna=substr($ref_dna, 1);
                                                        $alt_dna=substr($alt_dna, 1);
                                                }elsif(substr($ref_dna, 2, 1) eq substr($alt_dna, 2, 1))
                                                {
                                                        $ref_dna=substr($ref_dna, 0, 2);
                                                        $alt_dna=substr($alt_dna, 0, 2);
                                                }else
                                                {
                                                        die "something is wrong\n";
                                                }
                                        }
					$result=join("\t", @h)."\tMerged\=Yes\n".join("\t", @l[0..4])."\t$new_start\t$new_end\t$l[7]\t$effect\tMNP\t$ref_dna\t$alt_dna\t$alt_dna\t".join("\t",@l[13..33])."\tc\.$p$ref_dna\>$alt_dna\tp\.$codon_abbr{$ref_aa}$poa$codon_abbr{$alt_aa}\tp\.$ref_aa$poa$alt_aa\t".join("\t", @l[37..$#l])."\n";
					return $result;
				}
			}
		}
	}
}

###############=> END <=##################
##########################################


##########################################
########## Processing the data ###########
##############=> START <=#################
my $maf_tmp=$dirlog."/tmp";
system("sort -t'	' -k 1,1 -k 16,16 -k 6,6n $maf | grep -v ^# | grep -v ^Hugo > $maf_tmp");

#sort by gene, sample and start position
open(MAF, "$maf_tmp");

my $gene_sample = "";
my $candidates = "";
my $raw_output = "";
my ($start, $end);
while(<MAF>)
{
	chomp;
	my @line = split(/\t/, );
	if($snvonlyFlag)
	{
		next if $line[9] ne "SNP";
	}
	if($line[0].$line[15] ne $gene_sample)
	{
		if($candidates=~m/\n/)
		{
			$raw_output = readsConfirmed($candidates);
			if($raw_output)
			{
				print OUT "$raw_output";
			}
		}
		$candidates = $_;
	}else
	{
		if($line[5]-$end <= $distance)
		{
			$candidates.="\n$_";
		}else
		{
			if($candidates=~m/\n/)
			{
				$raw_output = readsConfirmed($candidates);
				if($raw_output)
				{
					print OUT "$raw_output";
				}
			}
			$candidates = $_;
		}
	}
	$gene_sample = $line[0].$line[15];
	$start = $line[5];
	$end = $line[6];   
}
if($candidates=~m/\n/)
{
	$raw_output = readsConfirmed($candidates);
	if($raw_output)
	{
		print OUT "$raw_output";
	}
}
close MAF;
close OUT;

if($mergeFlag)
{
	if(!$snvonlyFlag)
	{
		die "\"\-\-snvonly\" is required for MERGING mode\n";
	}else
	{
#		print STDERR "Warning: variants with distance > 2 will not be considered since they don't affect the same codon\n";
#		print STDERR "Merging...\n";
		open(MNP, "$output");
		my $header=<MNP>;
		print OUTPUT $header;
		my $datain = "";
		while(<MNP>)
		{
			chomp;
			if($_=~m/^#/)
			{
				if($datain ne "")
				{
					print OUTPUT toMerge($datain);
				}
				$datain=$_."\n";
			}else
			{
				$datain.=$_."\n";
			}
		}
		print OUTPUT toMerge($datain);
	}
}

## remove tmp ##
`rm $maf_tmp`;
`rm $f_cds_bed`;
`rm $f_gtf_cds`;
###############=> END <=##################
###########################################

