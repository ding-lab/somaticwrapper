sub bsub_pindel{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $pindel_dir = shift;
    my $sw_dir = shift;
    my $f_centromere = shift;

    $current_job_file = "j5_pindel".$sample_name.".sh";  
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $pindel_out = "$sample_full_path/pindel/pindel_out";
    system("mkdir -p $pindel_out");

    my $config_fn = "$pindel_out/$sample_name.config";
    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$IN_bam_T\t500\t$sample_name.T
$IN_bam_N\t500\t$sample_name.N
EOF

    my $pindel_args = " -T 4 -m 6 -w 1 ";

    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

$pindel_dir/pindel -f $REF -i $config_fn -o $pindel_out/$sample_name $pindel_args -J $f_centromere

EOF

    close OUT;

    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com );   
}

1;
