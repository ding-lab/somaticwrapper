```
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
```

# Email from Cyriac Kandoth 8/19/18:

By design, vcf2maf will never remove variants. The most common reason to have
fewer mutations in the output MAF is because the input had duplicate entries
for the same event+sample, or because an event was formatted incorrectly and
VEP fails on it. De-duplication happens silently, but the formatting issues are
reported as warnings/errors.

The "--filter-vcf" parameter currently only supports the ExAC VCF format. It
extracts ExAC allele counts from the VCF, and uses it to tag events as
"common_variant" in the FILTER column, as described here -
https://github.com/mskcc/vcf2maf/blob/v1.6.16/docs/vep_maf_readme.txt#L134-L138

I opted not to use the ExAC allele counts reported by VEP because for these
reasons -
https://github.com/mskcc/vcf2maf/blob/v1.6.16/data/known_somatic_sites.bed#L6-L9
- they include many known somatic hotspots that we never want to tag as
"common_variant". Instead, I construct my own variant of the
ExAC_nonTCGA.r0.3.1 VCF as described here -
https://gist.github.com/ckandoth/f265ea7c59a880e28b1e533a6e935697

Extra columns for gnomAD AFs were added to the MAF as of vcf2maf v1.6.16. But
it does not yet do any filtering/tagging using gnomAD data, because I haven't
had a chance to curate it for problems.
