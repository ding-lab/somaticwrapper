# note that f_exac is specific to GRCh37
# This is incomplete:
# * need to specify vep paths
# * need to figure out f_exac with GRCh38 issues

sub vcf_2_maf{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $f_exac = shift;

    my $bsub = "bash";
    $current_job_file = "j9_vcf_2_maf.".$sample_name.".sh";

    my $merged_results = "$sample_full_path/merged";
    my $filter_results = "$sample_full_path/maf";

    my $F_VCF_1="$merged_results/merged.vcf";
    my $F_VCF_2="$filter_results/$sample_name.vcf";
#    my $F_VEP_1="$merged_results/merged.VEP.vcf";
#    my $F_VEP_2="$filter_results/$sample_name.vep.vcf";

    my $F_maf="$filter_results/$sample_name.maf";
    system("ln -s $F_VCF_1 $F_VCF_2");
#    system("ln -s $F_VEP_1 $F_VEP_2");

    if (defined $f_exac) {
        my $exac_filter="--filter-vcf $f_exac";
    } else {
        my $exac_filter="";
    }

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
$perl $sw_dir/vcf2maf.pl --input-vcf $F_VCF_2 --output-maf $F_maf --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $REF $exac_filter

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 

}

1;
