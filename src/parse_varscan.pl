# The following files were created in $sample_full_path/varscan_out
#  bamfilelist.inp
#  varscan.out.som.log
#  varscan.out.som_indel.vcf 
#  varscan.out.som_snv.vcf 
# processing which takes place here will be written to varscan/filter_out 
#
# Principal output and CWL mapping:
    # varscan.out.som_snv.gvip.Somatic.hc.vcf      -> varscan_snv_process
    # varscan.out.som_indel.gvip.Somatic.hc.vcf    -> varscan_indel_process
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf   -> varscan_snv_filtered
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf   -> varscan_snv_dbsnp
    # varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf    -> varscan_indel_dbsnp

# CWL changes:
# * get rid of genomevip_label steps
# * input filenames are passed explicitly: varscan_indel_raw and varscan_snv_raw
#   * these had been varscan.out.som_indel.vcf and varscan.out.som_snv.vcf


sub parse_varscan{
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $dbsnp_db = shift;
    my $snpsift_jar = shift;
    my $varscan_jar = shift;
    my $varscan_indel_raw = shift; # indeloutgvip
    my $varscan_snv_raw = shift;  # snvoutgvip

    $current_job_file = "j4_parse_varscan.sh";

    my $bsub = "bash";
    my $filter_results = "$sample_full_path/varscan/filter_out";
    system("mkdir -p $filter_results");


    # These based on original script
    my $snvoutbase="$filter_results/varscan.out.som_snv";
    my $indeloutbase="$filter_results/varscan.out.som_indel";

    #my $thissnvorig="${snvoutbase}.gvip.Somatic.hc.vcf";  # This is genrated by varscan processSomatic
    my $thissnvorig="$filter_results/varscan.out.som_snv.gvip.Somatic.hc.vcf";  # This is genrated by varscan processSomatic
    my $myindelorig="${indeloutbase}.gvip.vcf";
    my $thissnvpass="${snvoutbase}.gvip.Somatic.hc.somfilter_pass.vcf";


    my $log_file="$filter_results/varscan.out.som.log";

    my $somsnvpass="$filter_results/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf";

    my $indeloutbase="$filter_results/varscan.out.som_indel";

    my $out = "$filter_results/vs_dbsnp_filter.snv.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.snv.annotator = $snpsift_jar
varscan.dbsnp.snv.db = $dbsnp_db
varscan.dbsnp.snv.rawvcf = $somsnvpass
varscan.dbsnp.snv.mode = filter
varscan.dbsnp.snv.passfile  = $filter_results/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf
varscan.dbsnp.snv.dbsnpfile = $filter_results/varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf
EOF

    my $out = "$filter_results/vs_dbsnp_filter.indel.input";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
varscan.dbsnp.indel.annotator = $snpsift_jar
varscan.dbsnp.indel.db = $dbsnp_db
varscan.dbsnp.indel.rawvcf = $indeloutbase.gvip.Somatic.hc.vcf
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
export JAVA_OPTS=\"-Xms256m -Xmx10g\"

echo \'APPLYING PROCESS FILTER TO SOMATIC SNVS:\' # &> $log_file
# Script below creates:
    # varscan.out.som_snv.gvip.Somatic.hc.vcf      -> used for SNV SNP filter below and vep annotation
    # varscan.out.som_snv.gvip.Somatic.vcf        
    # varscan.out.som_snv.gvip.LOH.hc.vcf         
    # varscan.out.som_snv.gvip.LOH.vcf            
    # varscan.out.som_snv.gvip.Germline.hc.vcf    
    # varscan.out.som_snv.gvip.Germline.vcf       
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $varscan_snv_raw $somatic_snv_params # &>> $log_file

echo \'APPLYING PROCESS FILTER TO SOMATIC INDELS:\' # &>> $log_file
# Script below creates:
    # varscan.out.som_indel.gvip.Germline.hc.vcf    
    # varscan.out.som_indel.gvip.Germline.vcf       
    # varscan.out.som_indel.gvip.LOH.hc.vcf         
    # varscan.out.som_indel.gvip.LOH.vcf            
    # varscan.out.som_indel.gvip.Somatic.hc.vcf     -> used for Indel SnP Filter below and vep annotation
    # varscan.out.som_indel.gvip.Somatic.vcf        
java \${JAVA_OPTS} -jar $varscan_jar processSomatic $varscan_indel_raw   $somatic_indel_params  # &>> $log_file


### Somatic Filter filters SNV based on indel
# http://varscan.sourceforge.net/using-varscan.html#v2.3_somaticFilter
echo \'APPLYING SOMATIC FILTER:\' # &>> $log_file

# Script below creates:
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf   -> used for SNV dbSnP and vep annotation
java \${JAVA_OPTS} -jar $varscan_jar somaticFilter  $thissnvorig $somatic_filter_params  --indel-file  $varscan_indel_raw --output-file  $somsnvpass  # &>> $log_file   

### dbSnP Filter

# 1) SNV
# Script below reads:  
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf
# and generates:
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_present.vcf  
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf     -> used for merge_vcf
    # varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_anno.vcf   
$perl $gvip_dir/dbsnp_filter.pl  $filter_results/vs_dbsnp_filter.snv.input

# 2) indel
# Script below reads
    # varscan.out.som_indel.gvip.Somatic.hc.vcf
# and generates:
    # varscan.out.som_indel.gvip.Somatic.hc.dbsnp_present.vcf 
    # varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf    -> used for merge_vcf
    # varscan.out.som_indel.gvip.Somatic.hc.dbsnp_anno.vcf    
$perl $gvip_dir/dbsnp_filter.pl $filter_results/vs_dbsnp_filter.indel.input

# Genome VIP SNV filter , and two GenomeVIP VEP annotation calls deleted from end of workflow per 
# discussion with Song

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 
}

1;
