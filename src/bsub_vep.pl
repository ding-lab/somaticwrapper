# By default, VEP searches for caches in $HOME/.vep; to use a different directory when running VEP, use --dir_cache.

# Dir cache: /data/D_VEP


my $assembly="GRCh37";
my $cachedir="/data/D_VEP";
#my $reffasta="$cachedir/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
my $reffasta="/data/A_Reference/demo20.fa";

# Note that VEP that had been used was v81, which used the script variant_effect_predictor.pl
# Newer versions of VEP use vep.pl as the command.  It remains to be seen what differences there are
# old: VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
my $vep_cmd="/usr/local/ensembl-vep/vep";

sub write_vep_input {
    my $config_fn = shift;
    my $module = shift;  # e.g. varscan.vep or strelka.vep
    my $vcf = shift;
    my $output = shift; 
    my $vep_cmd = shift;
    my $cache_dir = shift;
    my $reffasta = shift;
    my $assembly = shift;

    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$module.vcf = $vcf
$module.output = $output
$module.vep_cmd = $vep_cmd
$module.cachedir = $cachedir
$module.reffasta = $reffasta
$module.assembly = $assembly
EOF

}

sub bsub_vep{
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

    my $config_fn = "$filter_results/vs_vep.snv.input";
    my $vcf = "$varscan_results/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf";
    my $output = "$filter_results/varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf";
    my $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);


    $config_fn = "$filter_results/vs_vep.indel.input";
    $vcf = "$varscan_results/varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf";
    $output = "$filter_results/varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf";
    $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);

    $config_fn = "$filter_results/strelka_vep.snv.input";
    $vcf = "$strelka_results/strelka.somatic.snv.all.gvip.dbsnp_pass.vcf";
    $output = "$filter_results/strelka.somatic_snv.current_final.gvip.Somatic.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);

    $config_fn = "$filter_results/strelka_vep.indel.input";
    $vcf = "$strelka_results/strelka.somatic.indel.all.gvip.dbsnp_pass.vcf";
    $output = "$filter_results/strelka.somatic_indel.current_final.gvip.Somatic.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);

    $config_fn = "$filter_results/vs_vep.snv.inital.input";
    $vcf = "$varscan_results/varscan.out.som_snv.gvip.vcf";
    $output = "$filter_results/varscan.out.som_snv.gvip.VEP.vcf";
    $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);

    $config_fn = "$filter_results/vs_vep.indel.initial.input";
    $vcf = "$varscan_results/varscan.out.som_indel.gvip.vcf";
    $output = "$filter_results/varscan.out.som_indel.gvip.VEP.vcf";
    $module = "varscan.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);

    $config_fn = "$filter_results/strelka_vep.snv.initial.input";
    $vcf = "$strelka_results/strelka.somatic.snv.strlk_pass.gvip.vcf";
    $output = "$filter_results/strelka.somatic.snv.strlk_pass.gvip.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);

    $config_fn = "$filter_results/strelka_vep.indel.initial.input";
    $vcf = "$strelka_results/strelka.somatic.indel.strlk_pass.gvip.vcf";
    $output = "$filter_results/strelka.somatic.indel.strlk_pass.gvip.VEP.vcf";
    $module = "strelka.vep";
    write_vep_input($config_fn, $module, $vcf, $output, $vep_cmd, $cache_dir, $reffasta, $assembly);


    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

#export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8
#export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin
#export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64
export JAVA_OPTS=\"-Xms256m -Xmx512m\"
#export PATH=\${JAVA_HOME}/bin:\${PATH}

$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.snv.input >& $filter_results/vs_vep.snv.log
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.indel.input >& $filter_results/vs_vep.indel.log
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.snv.initial.input >& $filter_results/vs_vep.snv.initial.log
$perl $gvip_dir/vep_annotator.pl $filter_results/vs_vep.indel.initial.input >& $filter_results/vs_vep.indel.initial.log


$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.snv.input >& $filter_results/strelka_vep.snv.log
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.indel.input >& $filter_results/strelka_vep.indel.log
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.snv.initial.input >& $filter_results/strelka_vep.snv.initial.log
$perl $gvip_dir/vep_annotator.pl $filter_results/strelka_vep.indel.initial.input >& $filter_results/strelka_vep.indel.initial.log
EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    print("Aborting\n");
    die();

    system ( $bsub_com ); 


}

1;
