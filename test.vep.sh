#!/bin/bash
F_VCF_1=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.withmutect.vcf
F_VCF_1_filtered=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.filtered.withmutect.vcf
F_VCF_2=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/TCGA-2G-AAEQ-01.withmutect.vcf
F_VCF_2_filtered=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/TCGA-2G-AAEQ-01.withmutect.filtered.vcf
F_VEP_1=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.VEP.withmutect.vcf
F_VEP_1_filtered=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.VEP.withmutect.filtered.vcf
F_VEP_2=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/TCGA-2G-AAEQ-01.withmutect.vep.vcf
F_VEP_2_filtered=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/TCGA-2G-AAEQ-01.withmutect.filtered.vep.vcf
F_maf=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/TCGA-2G-AAEQ-01.withmutect.maf
F_maf_filtered=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/TCGA-2G-AAEQ-01.withmutect.filtered.maf
RUNDIR=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01
F_log=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/vep.merged.withmutect.log
merged_vep_vcf=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.withmutect.vcf
merged_vep_output=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.VEP.withmutect.vcf
rm ${F_log}
F_log_filtered=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/vep.merged.withmutect.filtered.log
#merged.vep.filtered.vcf=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.filtered.withmutect.vcf
#merged.vep.filtered.output=/gscmnt/gc2541/scao/tgct/somaticg36/TCGA-2G-AAEQ-01/merged.VEP.withmutect.filtered.vcf
rm ${F_log_filtered}
/gscmnt/gc2524/dinglab/scao/conda_root/envs/vep/bin/vep --species homo_sapiens --assembly GRCh38 --offline  --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir /gscmnt/gc2518/dinglab/scao/tools/vep/v102 --fasta /gscmnt/gc2518/dinglab/scao/tools/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa --input_file ${merged_vep_vcf} --output_file ${merged_vep_output} --polyphen b --af --af_1kg --af_esp --regulatory
