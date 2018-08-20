# Output is to pindel/pindel_out

# CWL addition:
# the `grep ChrID` step which was in the parse_pindel step (see SomaticWrapper.v2.Combined.pdf) has
# been moved to the run_pindel step
# The output of this step is then 
# * output port: pindel/pindel_out/pindel-raw.dat

# After running Pindel, pull out all reads from pindel raw output with the label ChrID
#   http://gmt.genome.wustl.edu/packages/pindel/user-manual.html

# runs pindel with arguments -T 4 -m 6 -w 1
# f_centromere is optional bed file passed to pindel to exclude regions
#   if it is defined, must exist
# will delete temporary files (pindel_D, etc) unless no_delete_temp is true

sub run_pindel {
    my $IN_bam_T = shift;
    my $IN_bam_N = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $pindel_dir = shift;
    my $f_centromere = shift;
    my $no_delete_temp = shift;

    if (! $no_delete_temp) {
        $no_delete_temp = 0; # avoid empty variables
    }

    my $centromere_arg = "";
    die "Error: Tumor BAM $IN_bam_T does not exist\n" if (! -e $IN_bam_T);
    die "Error: Tumor BAM $IN_bam_T is empty\n" if (! -s $IN_bam_T);
    die "Error: Normal BAM $IN_bam_N does not exist\n" if (! -e $IN_bam_N);
    die "Error: Normal BAM $IN_bam_N is empty\n" if (! -s $IN_bam_N);
    if ($f_centromere) {
        die "Error: Centromere BED file $f_centromere does not exist\n" if (! -e $f_centromere);
        $centromere_arg = "-J $f_centromere";
        print STDERR "Using centromere BED $f_centromere\n";
    } else {
        print STDERR "No centromere BED\n";
    }


    my $bsub = "bash";
    $current_job_file = "j5_pindel.sh";  
    my $step_output_fn = "pindel-raw.dat";

    my $pindel_out = "$sample_full_path/pindel/pindel_out";
    system("mkdir -p $pindel_out");

    my $config_fn = "$pindel_out/pindel.config";
    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$IN_bam_T\t500\tpindel.T
$IN_bam_N\t500\tpindel.N
EOF


    my $pindel_args = " -T 4 -m 6 -w 1 ";

# This step the original invocation from parse_pindel
#echo Collecting results in $pindel_results
#find $pindel_results -name \'*_D\' -o -name \'*_SI\' -o -name \'*_INV\' -o -name \'*_TD\'  > $outlist
#list=\$(xargs -a  $outlist)
#cat \$list | grep ChrID > $pin_var_file

    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

$pindel_dir/pindel -f $REF -i $config_fn -o $pindel_out/pindel $pindel_args $centromere_arg

cd $pindel_out 
grep ChrID pindel_D pindel_SI pindel_INV pindel_TD > $step_output_fn

if [[ $no_delete_temp == 1 ]]; then

>&2 echo Not deleting intermediate pindel files

else

>&2 echo Deleting intermediate pindel files
rm -f pindel_*

fi


EOF

    close OUT;

    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    my $return_code = system ( $bsub_com );
    die("Exiting ($return_code).\n") if $return_code != 0;
}

# NOTE from Sunantha, how to run per chrom in parallel

ORIG
echo "$IN_bam_T\t500\t$sample_name.T" > \${CONFIG}
echo "$IN_bam_N\t500\t$sample_name.N" >> \${CONFIG}
$pindel -T 4 -f $h37_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"." -m 6 -w 1 -J $f_centromere

NEW
echo "$IN_bam_T\t500\t$sample_name.T" > \${CONFIG}
echo "$IN_bam_N\t500\t$sample_name.N" >> \${CONFIG}
for chr in {1..22} X 
do 
nohup $pindel -T 4 -c \$chr -f $h37_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"."_\${chr}"." -m 6 -w 1 -J $f_centromere > \${myRUNDIR}"."/$sample_name"."_\${chr}_pindel.log"." & 
done 

PINDEL options:
           -c/--chromosome
           Which chr/fragment. Pindel will process reads for one chromosome each
           time. ChrName must be the same as in reference sequence and in read
           file. '-c ALL' will make Pindel loop over all chromosomes. The search
           for indels and SVs can also be limited to a specific region; -c
           20:10,000,000 will only look for indels and SVs after position
           10,000,000 = [10M, end], -c 20:5,000,000-15,000,000 will report
           indels in the range between and including the bases at position
           5,000,000 and 15,000,000 = [5M, 15M]. (default ALL)



1;
