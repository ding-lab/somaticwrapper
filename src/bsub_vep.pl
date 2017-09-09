
sub bsub_vep{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;

    $current_job_file = "j6_vep".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";


#cat > \${RUNDIR}/varscan/vs_vep.snv.input <<EOF
    my $out = "$filter_results/vs_vep.snv.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.vep.vcf = ./varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf
varscan.vep.output = ./varscan.out.som_snv.current_final.gvip.Somatic.VEP.vcf
varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
varscan.vep.assembly = GRCh37
EOF

#cat > \${RUNDIR}/varscan/vs_vep.indel.input <<EOF
    my $out = "$filter_results/vs_vep.indel.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.vep.vcf = ./varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf
varscan.vep.output = ./varscan.out.som_indel.current_final.gvip.Somatic.VEP.vcf
varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
varscan.vep.assembly = GRCh37
EOF

#cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.snv.input <<EOF
    my $out = "$filter_results/strelka_vep.snv.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
strelka.vep.vcf = ./strelka.somatic.snv.all.gvip.dbsnp_pass.vcf
strelka.vep.output = ./strelka.somatic_snv.current_final.gvip.Somatic.VEP.vcf
strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
strelka.vep.assembly = GRCh37
EOF

# cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.indel.input <<EOF
    my $out = "$filter_results/strelka_vep.indel.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
strelka.vep.vcf = ./strelka.somatic.indel.all.gvip.dbsnp_pass.vcf
strelka.vep.output = ./strelka.somatic_indel.current_final.gvip.Somatic.VEP.vcf
strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
strelka.vep.assembly = GRCh37
EOF

#cat > \${RUNDIR}/varscan/vs_vep.snv.inital.input <<EOF
    my $out = "$filter_results/vs_vep.snv.inital.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.vep.vcf = ./varscan.out.som_snv.gvip.vcf
varscan.vep.output = ./varscan.out.som_snv.gvip.VEP.vcf
varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
varscan.vep.assembly = GRCh37
EOF

#cat > \${RUNDIR}/varscan/vs_vep.indel.initial.input <<EOF
    my $out = "$filter_results/vs_vep.indel.initial.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.vep.vcf = ./varscan.out.som_indel.gvip.vcf
varscan.vep.output = ./varscan.out.som_indel.gvip.VEP.vcf
varscan.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
varscan.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
varscan.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
varscan.vep.assembly = GRCh37
EOF

#cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.snv.initial.input <<EOF
    my $out = "$filter_results/strelka_vep.snv.initial.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
strelka.vep.vcf = ./strelka.somatic.snv.strlk_pass.gvip.vcf
strelka.vep.output = ./strelka.somatic.snv.strlk_pass.gvip.VEP.vcf
strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
strelka.vep.assembly = GRCh37
EOF

#cat > \${RUNDIR}/strelka/strelka_out/results/strelka_vep.indel.initial.input <<EOF
    my $out = "$filter_results/strelka_vep.indel.initial.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
strelka.vep.vcf = ./strelka.somatic.indel.strlk_pass.gvip.vcf
strelka.vep.output = ./strelka.somatic.indel.strlk_pass.gvip.VEP.vcf
strelka.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl
strelka.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache
strelka.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
strelka.vep.assembly = GRCh37
EOF




    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash
scr_t0=\`date \+\%s\`
TBAM=".$sample_full_path."/".$sample_name.".T.bam
NBAM=".$sample_full_path."/".$sample_name.".N.bam
myRUNDIR=".$sample_full_path."/varscan
STATUSDIR=".$sample_full_path."/status
RUNDIR=".$sample_full_path."
export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8
export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin
export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64
export JAVA_OPTS=\"-Xms256m -Xmx512m\"
export PATH=\${JAVA_HOME}/bin:\${PATH}



cd \${RUNDIR}/varscan
$perl $gvip_dir/vep_annotator.pl ./vs_vep.snv.input >& ./vs_vep.snv.log
$perl $gvip_dir/vep_annotator.pl ./vs_vep.indel.input >& ./vs_vep.indel.log
$perl $gvip_dir/vep_annotator.pl ./vs_vep.snv.initial.input >& ./vs_vep.snv.initial.log
$perl $gvip_dir/vep_annotator.pl ./vs_vep.indel.initial.input >& ./vs_vep.indel.initial.log
cd \${RUNDIR}/strelka/strelka_out/results
$perl $gvip_dir/vep_annotator.pl ./strelka_vep.snv.input >& ./strelka_vep.snv.log
$perl $gvip_dir/vep_annotator.pl ./strelka_vep.indel.input >& ./strelka_vep.indel.log
$perl $gvip_dir/vep_annotator.pl ./strelka_vep.snv.initial.input >& ./strelka_vep.snv.initial.log
$perl $gvip_dir/vep_annotator.pl ./strelka_vep.indel.initial.input >& ./strelka_vep.indel.initial.log
EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    print("Aborting");
    die();

    system ( $bsub_com ); 


}

1;
