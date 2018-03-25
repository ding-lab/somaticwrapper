# CWL-specific changes
# * Get rid of unused arguments
# * pass all input VCFs
# * Output port: merged/merged.vcf

sub merge_vcf {
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $gatk_jar = shift;

    my $strelka_snv_vcf = shift;
    my $varscan_indel_vcf = shift;
    my $varscan_snv_vcf = shift;
    my $pindel_vcf = shift;

    my $bsub = "bash";
    $current_job_file = "j8_merge_vcf.sh";
    my $filter_results = "$sample_full_path/merged";
    system("mkdir -p $filter_results");


    my $merger_out = "$filter_results/merged.vcf";


    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xmx2g\"
 java \$JAVA_OPTS -jar $gatk_jar -R $REF -T CombineVariants -o $merger_out --variant:varscan $varscan_snv_vcf --variant:strelka $strelka_snv_vcf --variant:varindel $varscan_indel_vcf --variant:pindel $pindel_vcf -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,pindel,varindel


echo Written final result to $merger_out

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
