# Output is to pindel/pindel_out

sub run_pindel{
    my $IN_bam_T = shift;
    my $IN_bam_N = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $pindel_dir = shift;
    my $f_centromere = shift;

    my $bsub = "bash";
    $current_job_file = "j5_pindel.sh";  

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

    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

$pindel_dir/pindel -f $REF -i $config_fn -o $pindel_out/pindel $pindel_args -J $f_centromere

EOF

    close OUT;

    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com );   
}

1;
