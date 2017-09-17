my $assembly="GRCh37";
my $cachedir="/data/D_VEP";


sub merge_vcf {
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $gatk = shift;
    my $usedb = shift;  # 1 for testing/demo, 0 for production


    $current_job_file = "j8_merge_vcf.".$sample_name.".sh";
    my $filter_results = "$sample_full_path/merged";
    system("mkdir -p $filter_results");

    my $strelka_vcf = "$sample_full_path/strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf";
    my $varscan_vcf = "$sample_full_path/varscan/filter_out/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf";
    my $pindel_vcf = "$sample_full_path/pindel/filter_out/pindel.out.current_final.gvip.dbsnp_pass.vcf";
    my $varscan_indel = "$sample_full_path/varscan/filter_out/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf";
    my $merger_out = "$filter_results/merged.vcf";

#cat > \${RUNDIR}/vep.merged.input <<EOF
    my $out = "$filter_results/vep.merged.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
merged.vep.vcf = $merger_out
merged.vep.output = $filter_results/merged.VEP.vcf
merged.vep.vep_cmd = $vep_cmd
merged.vep.cachedir = $cachedir
merged.vep.reffasta = $REF
merged.vep.assembly = $assembly
merged.vep.usedb = $usedb
EOF

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xmx2g\"
# java \$JAVA_OPTS -jar $gatk -R $REF -T CombineVariants -o $merger_out --variant:varscan $varscan_vcf --variant:strelka $strelka_vcf --variant:varindel $varscan_indel --variant:pindel $pindel_vcf -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,pindel,varindel

$perl $gvip_dir/vep_annotator.pl $filter_results/vep.merged.input # &> $filter_results/vep.merged.log

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
