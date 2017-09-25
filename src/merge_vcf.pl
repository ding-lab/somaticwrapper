my $assembly="GRCh37";
my $cachedir="/data/D_VEP";

# NOTES on running 01BR001
# The file varscan/filter_out/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf has contig order which is 
# incompatible with the reference.  See https://software.broadinstitute.org/gatk/documentation/article.php?id=1328
# for how to correct this; to wit,
# java -jar picard.jar SortVcf I=original.vcf O=sorted.vcf SEQUENCE_DICTIONARY=reference.dict
# Where reference.dict is the .dict file associated with the reference.  Doing:
# java -jar /usr/local/picard.jar SortVcf I=varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf-orig O=varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf SEQUENCE_DICTIONARY=/data/A_Reference/hg19.dict
# Same problem with strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
# TODO: incorporate filtering of the following VCFs as the last filtering step:
# * strelka/filter_out/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
# * varscan/filter_out/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf
# * varscan/filter_out/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf
# * pindel/filter_out/pindel.out.current_final.gvip.dbsnp_pass.vcf


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
# db mode 1) uses online database (so cache isn't installed) 2) does not use tmp files
# It is meant to be used for testing and lightweight applications.  Use the cache for
# better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
    my $use_vep_db = shift;  # 1 for testing/demo, 0 for production


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
merged.vep.usedb = $use_vep_db
EOF

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xmx2g\"
 java \$JAVA_OPTS -jar $gatk -R $REF -T CombineVariants -o $merger_out --variant:varscan $varscan_vcf --variant:strelka $strelka_vcf --variant:varindel $varscan_indel --variant:pindel $pindel_vcf -genotypeMergeOptions PRIORITIZE -priority strelka,varscan,pindel,varindel

$perl $gvip_dir/vep_annotator.pl $filter_results/vep.merged.input  # &> $filter_results/vep.merged.log

echo Written final result to $merger_out

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
