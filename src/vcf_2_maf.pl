# note that f_exac is specific to GRCh37
# This is incomplete:
# * need to specify vep paths
# * need to figure out f_exac with GRCh38 issues

# vcf2maf.pl documentation 
# Usage:
#      perl vcf2maf.pl --help
#      perl vcf2maf.pl --input-vcf WD4086.vcf --output-maf WD4086.maf --tumor-id WD4086 --normal-id NB4086
# 
# Options:
#      --input-vcf      Path to input file in VCF format
#      --output-maf     Path to output MAF file
#      --tmp-dir        Folder to retain intermediate VCFs after runtime [Default: Folder containing input VCF]
#      --tumor-id       Tumor_Sample_Barcode to report in the MAF [TUMOR]
#      --normal-id      Matched_Norm_Sample_Barcode to report in the MAF [NORMAL]
#      --vcf-tumor-id   Tumor sample ID used in VCF's genotype columns [--tumor-id]
#      --vcf-normal-id  Matched normal ID used in VCF's genotype columns [--normal-id]
#      --custom-enst    List of custom ENST IDs that override canonical selection
#      --vep-path       Folder containing the vep script [~/vep]
#      --vep-data       VEP's base cache/plugin directory [~/.vep]
#      --vep-forks      Number of forked processes to use when running VEP [4]
#      --buffer-size    Number of variants VEP loads at a time; Reduce this for low memory systems [5000]
#      --any-allele     When reporting co-located variants, allow mismatched variant alleles too
#      --ref-fasta      Reference FASTA file [~/.vep/homo_sapiens/91_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz]
#      --filter-vcf     A VCF for FILTER tag common_variant. Set to 0 to disable [~/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz]
#      --max-filter-ac  Use tag common_variant if the filter-vcf reports a subpopulation AC higher than this [10]
#      --species        Ensembl-friendly name of species (e.g. mus_musculus for mouse) [homo_sapiens]
#      --ncbi-build     NCBI reference assembly of variants MAF (e.g. GRCm38 for mouse) [GRCh37]
#      --cache-version  Version of offline cache to use with VEP (e.g. 75, 84, 91) [Default: Installed version]
#      --maf-center     Variant calling center to report in MAF [.]
#      --retain-info    Comma-delimited names of INFO fields to retain as extra columns in MAF []
#      --min-hom-vaf    If GT undefined in VCF, minimum allele fraction to call a variant homozygous [0.7]
#      --remap-chain    Chain file to remap variants to a different assembly before running VEP
#      --help           Print a brief help message and quit
#      --man            Print the detailed manual

# TODO:
# Merge vcf_2_maf, annotate_vcf:
#   --Change --output_vep to --output_format, take on values [VCF, VEP, MAF]
#   --VCF and VEP processed with old annotate_vcf code, vcf_2_maf passed to Cyriac's code
# Make --cache-version, --assembly (= --ncbi-build) optional for MAF and VCF/VEP output
# Get rid of --cache_gz
# change --cache_dir to accept .tar.gz format, in which case it decompresses cache file in default location


# vars to set
--input-vcf
--output-maf
--vep-path
--vep-data - same as annotation?

sub vcf_2_maf{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $f_exac = shift;
    my $vcf2maf_dir = shift;
    my $merged_vcf = shift; # this is the output of merge step

    $current_job_file = "j9_vcf_2_maf.".$sample_name.".sh";

    my $merged_results = "$sample_full_path/merged";

    my $filter_results = "$sample_full_path/maf";
    system("mkdir -p $filter_results");

    my $out_maf = "$filter_results/merged.filtered.maf"; # this is the final output

    if (defined $f_exac) {
        my $exac_filter="--filter-vcf $f_exac";
    } else {
        my $exac_filter="";
    }

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
$perl /usr/local/mskcc-vcf2maf/vcf2maf.pl --input-vcf $merged_vcf --output-maf $out_maf --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $REF $exac_filter

EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com ); 

}

1;
