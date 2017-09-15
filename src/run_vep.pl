# By default, VEP searches for caches in $HOME/.vep; to use a different directory when running VEP, use --dir_cache.

# Dir cache: /data/D_VEP


# NOTE: the pipeline as implemented does not read teh output of this step, since annotaiton takes place
# after merge_vcf
# This code remains functional, however, and annotates various Varscan and Strelka (not pindel) output

my $assembly="GRCh37";
my $cachedir="/data/D_VEP";

sub write_vep_input {
    my $config_fn = shift;
    my $module = shift;  # e.g. varscan.vep or strelka.vep
    my $vcf = shift;
    my $output = shift; 
    my $vep_cmd = shift;
    my $cache_dir = shift;
    my $REF = shift;
    my $assembly = shift;
    my $vep_cmd = shift;

    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$module.vcf = $vcf
$module.output = $output
$module.vep_cmd = $vep_cmd
$module.cachedir = $cachedir
$module.reffasta = $REF
$module.assembly = $assembly
EOF

}

sub run_vep {
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $gvip_dir = shift;

    $current_job_file = "j6_vep".$sample_name.".sh";

    my $varscan_results = "$sample_full_path/varscan/filter_out";
    my $strelka_results = "$sample_full_path/strelka/filter_out";
    my $filter_results = "$sample_full_path/vep";
    system("mkdir -p $filter_results");

    # VCF exists, no data
    my $config_fn = "$filter_results/vs_vep.snv.input";
    my $vcf = "$varscan_results/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf";
    my $output = "$filter_results/varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf";
    my $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);

    # VCF exists, has data
    $config_fn = "$filter_results/vs_vep.indel.input";
    $vcf = "$varscan_results/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf";
    $output = "$filter_results/varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf";
    $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);

    # VCF exists, has data
    $config_fn = "$filter_results/vs_vep.snv.initial.input";   # formerly vs_vep.snv.inital.input, never ran
    $vcf = "$varscan_results/varscan.out.som_snv.gvip.vcf";
    $output = "$filter_results/varscan.out.som_snv.gvip.VEP.vcf";
    $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);

    # VCF exists, has data
    $config_fn = "$filter_results/vs_vep.indel.initial.input";
    $vcf = "$varscan_results/varscan.out.som_indel.gvip.vcf";
    $output = "$filter_results/varscan.out.som_indel.gvip.VEP.vcf";
    $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);

    # VCF exists, has data
    $config_fn = "$filter_results/strelka_vep.snv.input";
    $vcf = "$strelka_results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf";
    $output = "$filter_results/strelka.somatic_snv.current_final.gvip.Somatic.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);

    # VCF exists, has data
    $config_fn = "$filter_results/strelka_vep.indel.input";
    $vcf = "$strelka_results/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf";
    $output = "$filter_results/strelka.somatic_indel.current_final.gvip.Somatic.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);

    # VCF does not exist
    $config_fn = "$filter_results/strelka_vep.snv.initial.input";
    $vcf = "$strelka_results/strelka.somatic.snv.strlk_pass.gvip.vcf";
    $output = "$filter_results/strelka.somatic.snv.strlk_pass.gvip.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);

    # VCF does not exist 
    $config_fn = "$filter_results/strelka_vep.indel.initial.input";
    $vcf = "$strelka_results/strelka.somatic.indel.strlk_pass.gvip.vcf";
    $output = "$filter_results/strelka.somatic.indel.strlk_pass.gvip.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $REF, $assembly);


    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

export JAVA_OPTS=\"-Xms256m -Xmx512m\"

$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.snv.input &> $filter_results/vep.log
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.indel.input &>> $filter_results/vep.log
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.snv.initial.input &>> $filter_results/vep.log
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.indel.initial.input &>> $filter_results/vep.log


$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.snv.input &>> $filter_results/vep.log
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.indel.input &>> $filter_results/vep.log
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.snv.initial.input &>> $filter_results/vep.log
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.indel.initial.input &>> $filter_results/vep.log

echo Logs written to $filter_results/vep.log
EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
