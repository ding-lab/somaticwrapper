function index_vcf {
rm -f $1.gz $1.gz.tbi
bgzip $1
tabix -p vcf $1.gz
echo Written $1
}

mkdir -p dat

ORIGD=/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/somatic/test.for.matt/01BR001
NEWD=/gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data/data/01BR001

function run_strelka_evaluate {
OUT=$1
VCF=$2
cp $ORIGD/strelka/strelka_out/results/$VCF $OUT.orig.vcf
index_vcf $OUT.orig.vcf
 
cp $NEWD/strelka/strelka_out/results/$VCF $OUT.new.vcf
index_vcf $OUT.new.vcf
}

function parse_strelka_evaluate {
OUT=$1
VCF=$2
cp $ORIGD/strelka/strelka_out/results/$VCF $OUT.orig.vcf
index_vcf $OUT.orig.vcf

cp $NEWD/strelka/filter_out/$VCF $OUT.new.vcf
index_vcf $OUT.new.vcf
}

function run_varscan_evaluate {
OUT=$1
VCF=$2
cp $ORIGD/varscan/$VCF $OUT.orig.vcf
index_vcf $OUT.orig.vcf

cp $NEWD/varscan/varscan_out/$VCF $OUT.new.vcf
index_vcf $OUT.new.vcf
}

function parse_varscan_evaluate {
OUT=$1
VCF=$2
cp $ORIGD/varscan/$VCF $OUT.orig.vcf
index_vcf $OUT.orig.vcf

cp $NEWD/varscan/filter_out/$VCF $OUT.new.vcf
index_vcf $OUT.new.vcf
}

# 1) strelka run
#run_strelka_evaluate "dat/run_strelka.1a" "passed.somatic.snvs.vcf"
#run_strelka_evaluate "dat/run_strelka.1b" "passed.somatic.indels.vcf"

# 2) strelka filter
#parse_strelka_evaluate "dat/parse_strelka.2" "strelka.somatic.snv.all.gvip.dbsnp_pass.vcf"

# 3) varscan run
#run_varscan_evaluate "dat/run_varscan.3a" "varscan.out.som_snv.vcf"
#run_varscan_evaluate "dat/run_varscan.3b" "varscan.out.som_indel.vcf"

# 4) varscan filter
# steps 4.1, 4.2, 4.3 identical
parse_varscan_evaluate "dat/parse_varscan.4.1a" "varscan.out.som_snv.gvip.Somatic.hc.vcf"
parse_varscan_evaluate "dat/parse_varscan.4.1b" "varscan.out.som_snv.gvip.LOH.hc.vcf"
parse_varscan_evaluate "dat/parse_varscan.4.1c" "varscan.out.som_snv.gvip.Germline.hc.vcf"

parse_varscan_evaluate "dat/parse_varscan.4.2a" "varscan.out.som_indel.gvip.Germline.hc.vcf"
parse_varscan_evaluate "dat/parse_varscan.4.2b" "varscan.out.som_indel.gvip.LOH.hc.vcf"
parse_varscan_evaluate "dat/parse_varscan.4.2c" "varscan.out.som_indel.gvip.Somatic.hc.vcf"

parse_varscan_evaluate "dat/parse_varscan.4.3" "varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.vcf"

# Things get very different going forward
# TODO: confirm whether these are different/same compared to run with original dbsnp database
parse_varscan_evaluate "dat/parse_varscan.4.4" "varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
parse_varscan_evaluate "dat/parse_varscan.4.5" "varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf"


#parse_varscan_evaluate "dat/parse_varscan.4a" "varscan.out.som_snv.gvip.Somatic.hc.somfilter_pass.dbsnp_pass.vcf"
#parse_varscan_evaluate "dat/parse_varscan.4b" "varscan.out.som_indel.gvip.Somatic.hc.dbsnp_pass.vcf"
