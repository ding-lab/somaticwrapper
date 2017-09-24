my $assembly="GRCh37";
my $cachedir="/data/D_VEP";

my $snpsift_jar="/usr/local/snpEff/SnpSift.jar";

# filtered database created in B_Filter
my $datd="/data/B_Filter";
my $db="$datd/dbsnp.noCOSMIC.vcf.gz";
#my $db="$datd/short.dbsnp.noCOSMIC.vcf.gz";

# Skipping VEP annotation

# principal output used in merging: pindel.out.current_final.gvip.dbsnp_pass.vcf
# This is based on pindel.out.current_final.gvip.Somatic.vcf - this has just chr1
#   - who creates pindel.out.current_final.gvip.Somatic.vcf?

sub parse_pindel {
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $vep_cmd = shift;
    my $pindel_dir = shift;

    $current_job_file = "j7_parse_pindel".$sample_name.".sh";

    my $pindel_results = "$sample_full_path/pindel/pindel_out";
    my $filter_results = "$sample_full_path/pindel/filter_out";
    system("mkdir -p $filter_results");

    my $outlist="$filter_results/pindel.out.filelist";
    my $pin_var_file="$filter_results/pindel.out.raw";
    my $pre_current_final="$pin_var_file.CvgVafStrand_pass.Homopolymer_pass.vcf";
    my $current_final="$filter_results/pindel.out.current_final.gvip.Somatic.vcf";



## Pindel Filter
#cat > \${RUNDIR}/pindel/pindel_filter.input <<EOF
    my $out = "$filter_results/pindel_filter.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
pindel.filter.pindel2vcf = $pindel_dir/pindel2vcf
pindel.filter.variants_file = $pin_var_file
pindel.filter.REF = $REF
pindel.filter.date = 000000
pindel.filter.heterozyg_min_var_allele_freq = 0.2
pindel.filter.homozyg_min_var_allele_freq = 0.8
pindel.filter.mode = somatic
pindel.filter.apply_filter = true
pindel.filter.somatic.min_coverages = 10
pindel.filter.somatic.min_var_allele_freq = 0.05
pindel.filter.somatic.require_balanced_reads = \"true\"
pindel.filter.somatic.remove_complex_indels = \"true\"
pindel.filter.somatic.max_num_homopolymer_repeat_units = 6
EOF

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

# Note that in subsequent filtering (merge_vcf) only the file
#   pindel.out.current_final.gvip.dbsnp_pass.vcf
# is used, and the vep annotation is ignored.  
# For this reason, we're disabling VEP annotation here.

# #cat > \${RUNDIR}/pindel/pindel_vep.input <<EOF
#     my $module = "pindel.vep";
#     my $vcf = "$filter_results/pindel.out.current_final.gvip.dbsnp_pass.vcf";
#     my $output = "$filter_results/pindel.out.current_final.gvip.dbsnp_pass.VEP.vcf";
# 
#     my $out = "$filter_results/pindel_vep.input";
#     print("Writing to $out\n");
#     open(OUT, ">$out") or die $!;
#     print OUT <<"EOF";
# $module.vcf = $vcf
# $module.output = $output
# $module.vep_cmd = $vep_cmd
# $module.cachedir = $cachedir
# $module.reffasta = $REF
# $module.assembly = $assembly
# EOF

# 0. Pull out all reads from pindel raw output with the label ChrID
#   http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
# 1. run pindel_filter.  This produces
#    CvgVafStrand_pass 
#    CvgVafStrand_fail
#    Homopolymer_pass  
#    Homopolymer_fail  
# 2. Label things
# 3. Run dbSnP filter
# 4. Run VEP annotation  -> but we're skipping this because VEP annotation takes place downstream anyway

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash

find $pindel_results -name \'*_D\' -o -name \'*_SI\' -o -name \'*_INV\' -o -name \'*_TD\'  > $outlist
list=\$(xargs -a  $outlist)
cat \$list | grep ChrID > $pin_var_file
$perl $gvip_dir/pindel_filter.v0.5.pl $filter_results/pindel_filter.input

$perl $gvip_dir/genomevip_label.pl Pindel $pin_var_file.CvgVafStrand_pass.vcf $pin_var_file.CvgVafStrand_pass.gvip.vcf
$perl $gvip_dir/genomevip_label.pl Pindel $pre_current_final $pin_var_file.CvgVafStrand_pass.Homopolymer_pass.gvip.vcf 
$perl $gvip_dir/genomevip_label.pl Pindel $pin_var_file.CvgVafStrand_pass.Homopolymer_fail.vcf $pin_var_file.CvgVafStrand_pass.Homopolymer_fail.gvip.vcf 

cat $pin_var_file.CvgVafStrand_pass.Homopolymer_pass.gvip.vcf > $current_final

export JAVA_OPTS=\"-Xms256m -Xmx512m\"

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
