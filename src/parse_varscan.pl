my $datd="/data/B_Filter";
my $varscan_jar="/usr/local/VarScan.jar";
my $snpsift_jar="/usr/local/snpEff/SnpSift.jar";

# filtered database created in B_Filter
my $db="$datd/dbsnp.noCOSMIC.vcf.gz";
#my $db="$datd/short.dbsnp.noCOSMIC.vcf.gz";

# The following files were created in $sample_full_path/varscan
#  bamfilelist.inp
#  status  -  Note that this is empty, can probably be removed
#  varscan.out.som.log
#  varscan.out.som_indel.vcf
#  varscan.out.som_snv.vcf
# processing which takes place here will be written to $sample_full_path/varscan/filter_out ($filter_results)

sub parse_varscan{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;

    $current_job_file = "j4_parse_varscan".$sample_name.".sh";

    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $varscan_results = "$sample_full_path/varscan/varscan_out";
    my $filter_results = "$sample_full_path/varscan/filter_out";
    system("mkdir -p $filter_results");

    my $log_file="$filter_results/varscan.out.som.log";

    my $snvoutbase="$filter_results/varscan.out.som_snv";
    my $snvoutgvip="$snvoutbase.gvip.vcf";

    my $indeloutbase="$filter_results/varscan.out.som_indel";
    my $indeloutgvip="$indeloutbase.gvip.vcf";

    my $thissnvorig="$indeloutbase.gvip.Somatic.hc.vcf";
    my $thissnvpass="$indeloutbase.gvip.Somatic.hc.somfilter_pass.vcf";


# cat > $sample_full_path/varscan/vs_dbsnp_filter.snv.input <<EOF
    my $out = "$filter_results/vs_dbsnp_filter.snv.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.snv.annotator = $snpsift_jar
varscan.dbsnp.snv.db = $db
varscan.dbsnp.snv.rawvcf = $snvoutbase.gvip.Somatic.hc.somfilter_pass.vcf
#    should be according to Song's script: varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf
varscan.dbsnp.snv.mode = filter
varscan.dbsnp.snv.passfile  = $snvoutbase.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf
varscan.dbsnp.snv.dbsnpfile = $snvoutbase.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf
EOF

# cat > $sample_full_path/varscan/vs_dbsnp_filter.indel.input <<EOF
    my $out = "$filter_results/vs_dbsnp_filter.indel.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.indel.annotator = $snpsift_jar
varscan.dbsnp.indel.db = $db
varscan.dbsnp.indel.rawvcf = $thissnvorig
varscan.dbsnp.indel.mode = filter
varscan.dbsnp.indel.passfile  = $indeloutbase.gvip.Somatic.hc.dbsnp_pass.vcf
varscan.dbsnp.indel.dbsnpfile = $indeloutbase.gvip.Somatic.hc.dbsnp_present.vcf
EOF

    my $somatic_snv_params="--min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07";
    my $somatic_indel_params="--min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07";
    my $somatic_filter_params="--min-coverage 30 --min-reads2 4 --min-strands2 1 --min-avg-qual 20 --min-var-freq 0.10 --p-value 0.05";

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export VARSCAN_DIR="/usr/local"
# export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin
export JAVA_OPTS=\"-Xms256m -Xmx512m\"

# NOTE: comments below are to help track demo run.  Can be ignored and deleted.

# Script below creates varscan.out.som_snv.gvip.vcf
$perl $gvip_dir/genomevip_label.pl VarScan $varscan_results/varscan.out.som_snv.vcf $snvoutgvip
# Script below creates varscan.out.som_indel.gvip.vcf
$perl $gvip_dir/genomevip_label.pl VarScan $varscan_results/varscan.out.som_indel.vcf $indeloutgvip

echo \'APPLYING PROCESS FILTER TO SOMATIC SNVS:\' &> $log_file
# Script below creates files below in regular and "hc" versions
    # varscan.out.som_snv.gvip.Somatic.vcf  - not empty
    # varscan.out.som_snv.gvip.Germline.vcf - empty
    # varscan.out.som_snv.gvip.LOH.vcf      - not empty
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $snvoutgvip $somatic_snv_params &>> $log_file

echo \'APPLYING PROCESS FILTER TO SOMATIC INDELS:\' &>> $log_file
# Script below creates files below in regular and "hc" versions
    # varscan.out.som_indel.gvip.Somatic.vcf    - not empty
    # varscan.out.som_indel.gvip.Germline.vcf   - empty
    # varscan.out.som_indel.gvip.LOH.vcf        - empty
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $indeloutgvip   $somatic_indel_params  &>> $log_file

echo \'APPLYING SOMATIC FILTER:\' &>> $log_file
# Script below creates file,
    # varscan.out.som_indel.gvip.Somatic.hc.somfilter_pass.vcf  - empty
java \${JAVA_OPTS} -jar $varscan_jar somaticFilter  $thissnvorig $somatic_filter_params  --indel-file  $indeloutgvip --output-file  $thissnvpass  &>> $log_file

# Script below reads:
    # varscan.out.som_indel.gvip.Somatic.hc.somfilter_pass.vcf
# and generates:
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf  - empty
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf     - empty
    # varscan.out.som_indel.gvip.Somatic.hc.somfilter_pass.dbsnp_anno.vcf   - empty
$perl $gvip_dir/dbsnp_filter.pl  $filter_results/vs_dbsnp_filter.snv.input

# Script below reads
    # varscan.out.som_indel.gvip.Somatic.hc.vcf
# and generates:
    # varscan.out.som_indel.gvip.Somatic.hc.dbsnp_present.vcf - empty
    # varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf    - not empty
    # varscan.out.som_indel.gvip.Somatic.hc.dbsnp_anno.vcf    - not empty
$perl $gvip_dir/dbsnp_filter.pl $filter_results/vs_dbsnp_filter.indel.input

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
