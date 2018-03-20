

sub merge_vcf {
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $gatk = shift;
# db mode 1) uses online database (so cache isn't installed) 2) does not use tmp files
# It is meant to be used for testing and lightweight applications.  Use the cache for
# better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
    my $use_vep_db = shift;  # 1 for testing/demo, 0 for production
    my $output_vep = shift;  # output annotated vep rather than vcf format after merge step
    my $assembly = shift;
    my $cache_dir = shift;

    my $bsub = "bash";
    $current_job_file = "j8_merge_vcf.".$sample_name.".sh";
    my $filter_results = "$sample_full_path/merged";
    system("mkdir -p $filter_results");

    my $strelka_vcf = "$sample_full_path/strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf";
    my $varscan_vcf = "$sample_full_path/varscan/filter_out/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf";
    my $pindel_vcf = "$sample_full_path/pindel/filter_out/pindel.out.current_final.gvip.dbsnp_pass.vcf";
    my $varscan_indel = "$sample_full_path/varscan/filter_out/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf";
    my $merger_out = "$filter_results/merged.vcf";


    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xmx2g\"
 java \$JAVA_OPTS -jar $gatk -R $REF -T CombineVariants -o $merger_out --variant:varscan $varscan_vcf --variant:strelka $strelka_vcf --variant:varindel $varscan_indel --variant:pindel $pindel_vcf -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,pindel,varindel


echo Written final result to $merger_out

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
