
sub write_vep_input {
    my $config_fn = shift;
    my $module = shift;  # e.g. varscan.vep or strelka.vep
    my $vcf = shift;
    my $output = shift; 
    my $vep_cmd = shift;
    my $cache_dir = shift;   
    my $REF = shift;
    my $assembly = shift;
    my $use_vep_db = shift;  # 1 for testing/demo, 0 for production
    my $output_vep = shift;  # output annotated vep rather than vcf format after merge step.  add suffix 'vep' to output

    if ($output_vep) {
        $output = "$output.vep";
    }

    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$module.vcf = $vcf
$module.output = $output
$module.vep_cmd = $vep_cmd
$module.cachedir = $cache_dir
$module.reffasta = $REF
$module.assembly = $assembly
$module.usedb = $use_vep_db  
$module.output_vep = $output_vep  
EOF
}

sub run_vep {
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $assembly = shift;
    my $cache_dir = shift;
    my $use_vep_db = shift;  # 1 for testing/demo, 0 for production
    my $output_vep = shift;  # output annotated vep rather than vcf format after merge step
    my $annotate_intermediate = shift;  # annotate intermediate files (testing)

    $current_job_file = "j6_vep".$sample_name.".sh";

    my $varscan_results = "$sample_full_path/varscan/filter_out";
    my $strelka_results = "$sample_full_path/strelka/filter_out";
    my $pindel_results = "$sample_full_path/pindel/filter_out";
    my $merge_results = "$sample_full_path/merged";
    my $filter_results = "$sample_full_path/vep";
    system("mkdir -p $filter_results");

    if ($annotate_intermediate) {
        write_vep_input(
            "$filter_results/vs_vep.snv.input", 
            "varscan.vep",
            "$varscan_results/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf",
            "$filter_results/varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/vs_vep.indel.input",
            "varscan.vep",
            "$varscan_results/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf",
            "$filter_results/varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/vs_vep.snv.initial.input",
            "varscan.vep",
            "$varscan_results/varscan.out.som_snv.gvip.vcf",
            "$filter_results/varscan.out.som_snv.gvip.VEP.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/vs_vep.indel.initial.input",
            "varscan.vep",
            "$varscan_results/varscan.out.som_indel.gvip.vcf",
            "$filter_results/varscan.out.som_indel.gvip.VEP.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/strelka_vep.snv.input",
            "$strelka_results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf",
            "$filter_results/strelka.somatic_snv.current_final.gvip.Somatic.VEP.vcf",
            "strelka.vep",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

    ### This is not being calculated in parse_strelka 
    #    $config_fn = "$filter_results/strelka_vep.indel.input";
    #    $vcf = "$strelka_results/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf";  ###
    #    $output = "$filter_results/strelka.somatic_indel.current_final.gvip.Somatic.VEP.vcf";
    #    $module = "strelka.vep";
    #    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/strelka_vep.snv.initial.input",
            "strelka.vep",
            "$strelka_results/strelka.somatic.snv.strlk_pass.gvip.vcf",
            "$filter_results/strelka.somatic.snv.strlk_pass.gvip.VEP.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/strelka_vep.indel.initial.input",
            "strelka.vep",
            "$strelka_results/strelka.somatic.indel.strlk_pass.gvip.vcf",
            "$filter_results/strelka.somatic.indel.strlk_pass.gvip.VEP.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/pindel.initial.input",
            "pindel.vep",
            "$pindel_results/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.gvip.vcf",
            "$filter_results/pindel.initial.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

        write_vep_input(
            "$filter_results/pindel.final.input",
            "pindel.vep",
            "$pindel_results/pindel.out.current_final.gvip.dbsnp_pass.vcf",
            "$filter_results/pindel.final.vcf",
            $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);
    }

# Adding final output to vep annotation
    write_vep_input(
        "$filter_results/vep.merged.input",
        "merged.vep",
        "$merge_results/merged.vcf",
        "$filter_results/merged.VEP.vcf",
        $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx512m\"

$perl $gvip_dir/vep_annotator.pl $filter_results/vep.merged.input  

EOF

    if ($annotate_intermediate) {

        print OUT <<"EOF";
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.snv.input 
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.indel.input 
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.snv.initial.input 
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.indel.initial.input 

$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.snv.input 
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.indel.input 
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.snv.initial.input 
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.indel.initial.input 

$perl $gvip_dir/vep_annotator.pl $filter_results/pindel.initial.input 
$perl $gvip_dir/vep_annotator.pl $filter_results/pindel.final.input
EOF
    }

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
