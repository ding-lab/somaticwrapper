Set up Strelka demo.

1. Download test data distributed with Strelka
These are small BAMs which have variants in chrom "demo20", and an associated reference

2. Modify and stage the BAMs and reference for use in SomaticWarapper
* The unusual chrom name "demo20" cannot be annotated with VEP because it does not exist in
the database. For it to work, rename this as chrom "20"
* BAMs are renamed according to SomaticWrapper specifications

3. Create dbSnP filter
* The normal way of doing this is in B_filter, which downloads dbSnP and Cosmic databases
and filters them.  This is a long process because the datasets are large.  To make the demo
easier, we distribute a small version of the dbSnP and COSMIC database (created by B_Filter/4_makeStrelkaTestData.sh)
The filter processing is essentially the same as B_Filter/3_make_variant_filter.sh
