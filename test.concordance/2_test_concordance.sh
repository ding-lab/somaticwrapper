#TEST="bash ./snpSift_concordance.sh"
TEST="/gsc/bin/vcf-compare"

# vcf-compare -H A.vcf.gz B.vcf.gz C.vcf.gz

function run_test {
A=$1.orig.vcf.gz
B=$1.new.vcf.gz

echo Comparing $A and $B using $TEST

$TEST $A $B
}

#run_test dat/run_strelka.1a
#run_test dat/run_strelka.1b
#run_test dat/parse_strelka.2

#run_test dat/run_varscan.3a
#run_test dat/run_varscan.3b

run_test dat/parse_varscan.4.1a
run_test dat/parse_varscan.4.1b
run_test dat/parse_varscan.4.1c

run_test dat/parse_varscan.4.2a
run_test dat/parse_varscan.4.2b
run_test dat/parse_varscan.4.2c


run_test dat/parse_varscan.4.3
run_test dat/parse_varscan.4.4
run_test dat/parse_varscan.4.5
