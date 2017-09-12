# Test case
# running just vep
perl somatic_calling_v1.1.pl /data/data 6

# Issues with this run:
# 1. input data has zero variants, which creates problems with vep when header is stripped
#   Confirm why no variants in input file.  Expected?
# 2. strelka demo data has a fake chromosome "demo20".  This chrom is not found in the cache so it causes errors
#   TODO: Rename "demo20" -> "20" in all data
