
# Skipping VEP annotation

# principal output used in merging: pindel.out.current_final.gvip.dbsnp_pass.vcf
# This is based on pindel.out.current_final.gvip.Somatic.vcf 

sub parse_pindel {
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $pindel_dir = shift;
    my $db = shift;
    my $snpsift_jar = shift;
    my $pindel_config = shift;

    $current_job_file = "j7_parse_pindel.sh";

    my $bsub = "bash";
    my $pindel_results = "$sample_full_path/pindel/pindel_out";
    my $filter_results = "$sample_full_path/pindel/filter_out";
    system("mkdir -p $filter_results");

    my $outlist="$filter_results/pindel.out.filelist";
    my $pin_var_file="$filter_results/pindel.out.raw";
    #my $pre_current_final="$pin_var_file.CvgVafStrand_pass.Homopolymer_pass.vcf";
    my $pre_current_final="$filter_results/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf";
    my $current_final="$filter_results/pindel.out.current_final.gvip.Somatic.vcf";



## Pindel Filter - below is input into pindel_filter.v0.5
# lines below are added to data from $pindel_config
    die "$pindel_config does not exist\n" unless (-f $pindel_config);

    my $out = "$filter_results/pindel_filter.input";
    print("Copying $pindel_config to $out and appending\n");
    system("cp $pindel_config $out");

    open(OUT, ">>$out") or die $!;
    print OUT <<"EOF";
pindel.filter.pindel2vcf = $pindel_dir/pindel2vcf
pindel.filter.variants_file = $pin_var_file
pindel.filter.REF = $REF
pindel.filter.date = 000000
EOF

#pindel.filter.heterozyg_min_var_allele_freq = 0.2
#pindel.filter.homozyg_min_var_allele_freq = 0.8
#pindel.filter.mode = somatic
#pindel.filter.apply_filter = true
#pindel.filter.somatic.min_coverages = 10
#pindel.filter.somatic.min_var_allele_freq = 0.10
#pindel.filter.somatic.require_balanced_reads = \"true\"
#pindel.filter.somatic.remove_complex_indels = \"true\"
#pindel.filter.somatic.max_num_homopolymer_repeat_units = 6

## dbSnP Filter
#cat > \${RUNDIR}/pindel/pindel_dbsnp_filter.indel.input <<EOF
    my $out = "$filter_results/pindel_dbsnp_filter.indel.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
pindel.dbsnp.indel.annotator = $snpsift_jar
pindel.dbsnp.indel.db = $db
pindel.dbsnp.indel.rawvcf = $filter_results/pindel.out.current_final.gvip.Somatic.vcf
pindel.dbsnp.indel.mode = filter
pindel.dbsnp.indel.passfile  = $filter_results/pindel.out.current_final.gvip.dbsnp_pass.vcf
pindel.dbsnp.indel.dbsnpfile = $filter_results/pindel.out.current_final.gvip.dbsnp_present.vcf
EOF

# 0. Pull out all reads from pindel raw output with the label ChrID
#   http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
# 1. run pindel_filter.  This produces
#    pindel.out.raw.CvgVafStrand_pass 
#    pindel.out.raw.CvgVafStrand_fail
#    pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass  
#    pindel.out.raw.CvgVafStrand_pass.Homopolymer_fail  
# 2. Label things
# 3. Run dbSnP filter
# 4. Run VEP annotation  -> but we're skipping this because VEP annotation takes place downstream anyway

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash

echo Collecting results in $pindel_results
find $pindel_results -name \'*_D\' -o -name \'*_SI\' -o -name \'*_INV\' -o -name \'*_TD\'  > $outlist
list=\$(xargs -a  $outlist)
cat \$list | grep ChrID > $pin_var_file

echo Running pindel_filter.v0.5.pl
$perl $gvip_dir/pindel_filter.v0.5.pl $filter_results/pindel_filter.input

echo Running genomevip_label.pl
$perl $gvip_dir/genomevip_label.pl Pindel $filter_results/pindel.out.raw.CvgVafStrand_pass.vcf $filter_results/pindel.out.raw.CvgVafStrand_pass.gvip.vcf
#$perl $gvip_dir/genomevip_label.pl Pindel $pre_current_final $filter_results/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.gvip.vcf 
$perl $gvip_dir/genomevip_label.pl Pindel $filter_results/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf $filter_results/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.gvip.vcf 
$perl $gvip_dir/genomevip_label.pl Pindel $filter_results/pindel.out.raw.CvgVafStrand_pass.Homopolymer_fail.vcf $filter_results/pindel.out.raw.CvgVafStrand_pass.Homopolymer_fail.gvip.vcf 

# how does this differ from cp?
cat $pin_var_file.CvgVafStrand_pass.Homopolymer_pass.gvip.vcf > $current_final      # pindel.out.current_final.gvip.Somatic.vcf

export JAVA_OPTS=\"-Xms256m -Xmx10g\"

echo Running dbsnp_filter.pl
$perl $gvip_dir/dbsnp_filter.pl $filter_results/pindel_dbsnp_filter.indel.input

# Skipping VEP annotation because it is ignored in merge_vcf
# $perl $gvip_dir/vep_annotator.pl $filter_results/pindel_vep.input &> $filter_results/pindel_vep.log

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
