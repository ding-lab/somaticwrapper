# VAF Filter development

Goal: reproduce filtering by `vaf_filter_v1.2.pl` using python script
Part of effort to reproduce SomaticWrapper functionality with TinDaisy

## Development notes

Using vcf_filter.py from pyvcf seems to work:
```
vcf_filter.py --no-filtered varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf dps --depth-per-sample 100
```
It does not work on the merged filtered vcf for some reason.

## Understanding merging

Merge step, implemented in `somaticwrapper.cwl/src/merge_vcf.pl`, uses the following [CombineVariants] step,
```
java ... -T CombineVariants \
    --variant:varscan $varscan_snv_vcf \
    --variant:strelka $strelka_snv_vcf \
    --variant:varindel $varscan_indel_vcf \
    --variant:pindel $pindel_vcf \
    -genotypeMergeOptions PRIORITIZE \
    -priority strelka,varscan,pindel,varindel
```

Note, however, that the tumor, normal sample names in varscan and strelka are 'TUMOR', 'NORMAL', resp., while
for pindel they are `pindel.T` and `pindel.N`.  This keeps merge from combining the pindel and varindel calls.

Also, the merged pindel calls then take on the name of the sample:
```
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  C3N-01649.N C3N-01649.T NORMAL  TUMOR
```

To resolve this, it will be necessary to pre-process pindel results to make sure they have a standard naming convention
```
awk 'BEGIN{FS="\t";OFS="\t"}{if ($1 == "#CHROM") print $1, $2, $3, $4, $5, $6, $7, $8, $9, "NORMAL", "TUMOR"; else print}' pindel.short.vcf
```

[CombineVariants]: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php#



# Getting stated

Doing testing from within `cgc-images.sbgenomics.com/m_wyczalkowski/somatic-wrapper:cwl` container.  This has python
and other necessary packages installed.  Start container with,
```
bash start_docker.sh
```
This mounts the entire test project `/Users/mwyczalk/Projects/VCF_Process` as `/data` and allows docker work to make modifications
on the host.  Based on `/Users/mwyczalk/Projects/Rabix/TinDaisy/StrelkaDemo.docker.testing`.

```
git pull origin cwl  # Get latest image of SomaticWrapper CWL branch. Not necessary here 
cd /data
bash 1_run_VCF_Filter.py
```

Enter container from another terminal:
```
docker exec -it d0dec210bf3 bash
```

# Data

## Somatic Wrapper

Using `somaticwrapper.master/vaf_filter_v1.2.pl` as guide for implementation.

Most recent version of above cloned with,
`git clone https://github.com/ding-lab/somaticwrapper`
Subsequent edits made for clarity.

`vaf_filter_v1.2.pl (98c6a82)` has following filter logic, per Song Cao:

* minimum coverage filtering: 20
* tumor >= 5% and normal <=2% for varscan and strelka (percentages refer to VAF)
  * pindel tumor >=10% since the vaf calculation underestimates the ref coverage 
* filter for indel length (100 bps)  
* SNV: called by both strelka and varscan 
* INDEL: called by either varscan or pindel 

* 

## Somatic Wrapper CWL
```
git clone -b cwl https://github.com/ding-lab/somaticwrapper somaticwrapper.cwl
```

## VAF Calculations

VAF defined as (number of reads supporting variant)/(number of reads supporting variant and number of reads supporting reference) 
See also [here](https://www.biostars.org/p/226897/).




## Past Python work

https://github.com/ding-lab/MuSic_Remix/blob/master/calculate_mutation_rate.py
/Users/mwyczalk/Data/Virus/TCGA_Virus/BPS.TCGA_Virus.Lite/bps-core/src/util/processVCF.py  ** use this for basis of work
/Users/mwyczalk/Data/Virus/TCGA_Virus/BPS.TCGA_Virus.Lite/bps-core/src/analysis/vafFilter.py

# Data

## SomaticWrapper run

*Original runs (MGI)*
```
Results: /gscmnt/gc2533/dinglab/scao/cptac3/batch3/wxs/adjsomatic/C3N-01649
    top-level results (merge, vcf, etc) here: merged.C3N-01649.tar.gz
Run scripts: /gscmnt/gc2533/dinglab/scao/cptac3/batch3/wxs/tmpsomatic
    tmpsomatic.C3N-01649.tar.gz
Logs: /gscmnt/gc2533/dinglab/scao/cptac3/batch3/wxs/LSF_DIR_SOMATIC
    LSF_DIR_SOMATIC.C3N-01649.tar.gz
```

These are copied to `origdata` and expanded in `origdata/sw.C3N-01694`

*Results copy (Denali)*
Select SomaticWrapper results: denali:/diskmnt/Projects/Users/hsun/beta_tinDaisy/compare/mgi_sw_C3N-01649/

## Tin Daisy run (Denali)

Results: `denali:/diskmnt/Projects/Users/hsun/beta_tinDaisy/tin-daisy/results/TinDaisy.workflow-2018-07-30-135946.799/root/`


### Intermediates
```
scp mwyczalk_test@10.22.24.1:/diskmnt/Projects/Users/hsun/beta_tinDaisy/tin-daisy/results/TinDaisy.workflow-2018-07-30-135946.799/root/s3_parse_strelka/results/strelka/filter_out/strelka.somatic.snv.all.dbsnp_pass.vcf .
scp mwyczalk_test@10.22.24.1:/diskmnt/Projects/Users/hsun/beta_tinDaisy/tin-daisy/results/TinDaisy.workflow-2018-07-30-135946.799/root/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_snv.Somatic.hc.somfilter_pass.dbsnp_pass.vcf .
scp mwyczalk_test@10.22.24.1:/diskmnt/Projects/Users/hsun/beta_tinDaisy/tin-daisy/results/TinDaisy.workflow-2018-07-30-135946.799/root/s4_parse_varscan/results/varscan/filter_out/varscan.out.som_indel.Somatic.hc.dbsnp_pass.vcf .
scp mwyczalk_test@10.22.24.1:/diskmnt/Projects/Users/hsun/beta_tinDaisy/tin-daisy/results/TinDaisy.workflow-2018-07-30-135946.799/root/s7_parse_pindel/results/pindel/filter_out/pindel.out.current_final.dbsnp_pass.vcf .
```

### Final
`scp mwyczalk_test@10.22.24.1:/diskmnt/Projects/Users/hsun/beta_tinDaisy/tin-daisy/results/TinDaisy.workflow-2018-07-30-135946.799/root/s8_merge_vcf/results/merged/merged.vcf origdata/td.merged.vcf`

## Somatic Wrapper
```
scp mwyczalk_test@10.22.24.1:/diskmnt/Projects/Users/hsun/beta_tinDaisy/compare/mgi_sw_C3N-01649/merged.vcf origdata/sw.merged.vcf
scp mwyczalk_test@10.22.24.1:/diskmnt/Projects/Users/hsun/beta_tinDaisy/compare/mgi_sw_C3N-01649/merged.filtered.vcf origdata/sw.merged.filtered.vcf
```


# Analysis

## `set` values from sw.merged.vcf
```
grep -v "^#" sw.merged.vcf | cut -f 8 | tr ';' '\n' | grep "set=" | sort | uniq -c
 254 set=pindel
  24 set=strelka
   9 set=strelka-varscan
  88 set=varindel
1554 set=varscan
```

```
grep -v "^#" origdata/sw.merged.filtered.vcf | cut -f 8 | tr ';' '\n' | grep "set=" | sort | uniq -c
  15 set=pindel
   3 set=strelka-varscan
  40 set=varindel
```
