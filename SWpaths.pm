# parameters below based on Docker image locations.  It would perhaps be useful to define these in a configuration file.
package SWpaths;

$sw_dir = "/usr/local/somaticwrapper";
$gvip_dir="$sw_dir/src/GenomeVIP";
$filter_dir="$sw_dir/src/vcf_filters";
$strelka_dir = "/usr/local/strelka";
$strelka2_dir = "/usr/local/strelka2";
$gatk_jar = "/usr/local/GenomeAnalysisTK/GenomeAnalysisTK.jar";
$perl = "/usr/bin/perl";
$vep_cmd = "/usr/local/ensembl-vep/vep";
$pindel_dir = "/usr/local/pindel";
$snpsift_jar = "/usr/local/snpEff/SnpSift.jar";
$varscan_jar = "/usr/local/VarScan.jar";
$samtools="/usr/local/bin/samtools";
$vcf2maf="/usr/local/mskcc-vcf2maf/vcf2maf.pl";

1;
