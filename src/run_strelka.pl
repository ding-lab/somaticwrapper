# Run Strelka

# files written strelka/strelka_out/results: 
# * all.somatic.indels.vcf     
# * all.somatic.snvs.vcf
# * passed.somatic.indels.vcf   -> output
# * passed.somatic.snvs.vcf     -> output

sub run_strelka {
    my $IN_bam_T = shift;
    my $IN_bam_N = shift;
    my $results_dir = shift;
    my $job_files_dir = shift;
    my $STRELKA_DIR = shift;
    my $REF = shift;
    my $strelka_config = shift;

    my $bsub = "bash";
    $current_job_file = "j1_streka.sh"; 
    die "Error: Tumor BAM $IN_bam_T does not exist\n" if (! -e $IN_bam_T);
    die "Error: Tumor BAM $IN_bam_T is empty\n" if (! -s $IN_bam_T);
    die "Error: Normal BAM $IN_bam_N does not exist\n" if (! -e $IN_bam_N);
    die "Error: Normal BAM $IN_bam_N is empty\n" if (! -s $IN_bam_N);

#    my $strelka_config = "/usr/local/somaticwrapper/config/strelka.ini";
    my $strelka_out=$results_dir."/strelka/strelka_out";
    my $strelka_bin="$STRELKA_DIR/bin/configureStrelkaWorkflow.pl";


    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash

if [ -d $strelka_out ] ; then
    rm -rf $strelka_out
fi

$strelka_bin --normal $IN_bam_N --tumor $IN_bam_T --ref $REF --config $strelka_config --output-dir $strelka_out

cd $strelka_out
make -j 16
EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";

    print($bsub_com."\n");
    system ( $bsub_com );
}

1;
