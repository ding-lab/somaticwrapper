# Testing pyvcf's exensible vcf_flter.py framework
# looking at merged data.  Here, pindel samples != (strelka, varscan) samples

DATAD="/data"  # this works if running inside of docker

MERGED_VCF="$DATAD/dat/sw.merged.short.vcf"
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	C3N-01649.N	C3N-01649.T	NORMAL	TUMOR
#1	1229153	.	G	C	.	PASS	AC=1;AF=0.250;AN=4;DP=242;GENOMEVIP=1;GPV=1E0;SOMATIC;SPV=4.8701E-2;SS=2;SSC=13;set=varscan	GT:AD:DP:DP4:FREQ:RD	./.	./.	0/0:1:96:68,24,1,0:1.08%:92	0/1:9:146:87,49,8,1:6.21%:136
#1	14106395	.	C	CTCA	.	PASS	AC=0;AF=0.00;AN=4;END=14106395;GENOMEVIP=2;HOMLEN=2;HOMSEQ=TC;SVLEN=3;SVTYPE=INS;set=pindel	GT:AD	0/0:11,0	0/0:16,2	./.	./.
#1	20966723	.	CA	C	.	PASS	AC=0;AF=0.00;AN=4;END=20966724;GENOMEVIP=2;HOMLEN=4;HOMSEQ=AAAA;SVLEN=-1;SVTYPE=DEL;set=pindel	GT:AD	0/0:20,0	0/0:22,2	./.	./.
#1	24671747	.	GC	G	.	PASS	AC=1;AF=0.250;AN=4;DP=74;GENOMEVIP=1;GPV=1E0;SOMATIC;SPV=3.8913E-2;SS=2;SSC=14;set=varindel	GT:AD:DP:DP4:FREQ:RD	./.	./.	0/0:1:39:4,34,0,1:2.56%:38	0/1:6:35:1,28,0,6:17.14%:29
#1	31415038	.	TTATATATATATATATA	T	.	PASS	AC=1;AF=0.250;AN=4;DP=76;GENOMEVIP=1;GPV=1E0;SOMATIC;SPV=1.4021E-2;SS=2;SSC=18;set=varindel	GT:AD:DP:DP4:FREQ:RD	./.	./.	0/0:1:35:0,34,0,1:2.86%:34	0/1:9:41:0,32,0,9:21.95%:32
#1	145311533	.	G	C	.	PASS	AC=0;AF=0.00;AN=0;GENOMEVIP=5;NT=ref;QSS=16;QSS_NT=16;SGT=GG->CG;SOMATIC;TQSS=1;TQSS_NT=1;set=strelka	AU:CU:DP:FDP:GU:SDP:SUBDP:TU			0,0:0,1:25:0:25,83:0:0:0,0	0,0:6,9:45:0:38,109:0:0:1,1
#12	50745858	.	T	G	.	PASS	AC=0;AF=0.00;AN=0;DP=498;GPV=1E0;NT=ref;QSS=17;QSS_NT=17;SGT=TT->GT;SOMATIC;SPV=2.552E-2;SS=2;SSC=15;TQSS=2;TQSS_NT=2;set=strelka-varscan	AU:CU:DP:FDP:GU:SDP:SUBDP:TU			0,0:0,0:251:1:13,16:0:0:237,242	0,1:0,0:340:3:22,27:0:0:315,331

export PYTHONPATH="somaticwrapper.cwl/vcf_filters:$PYTHONPATH"
VAF_FILTER_LOCAL="vaf_filter.py"  # filter module

MAIN_FILTER="vcf_filter.py --no-filtered" # Assuming in path
CALLER="--caller merged"

# Arguments to VAF filter
SNV_VAF_ARGS="vaf --min_vaf_somatic 0.1 --debug" # --debug"
$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $MERGED_VCF $SNV_VAF_ARGS $CALLER


## testing - varscan SNP - OK
#CALLER="--caller varscan"
#$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $VARSCAN_VCF $SNV_VAF_ARGS $CALLER

# testing - varscan INDEL - OK
#CALLER="--caller varscan"
#$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $VARSCAN_INDEL_VCF $SNV_VAF_ARGS $CALLER

## testing - pindel.  Note that pindel has different sample names for now
#CALLER="--caller pindel" - OK
#NAMES="--normal_name pindel.N --tumor_name pindel.T"
#$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $PINDEL_VCF $SNV_VAF_ARGS $CALLER $NAMES

# testing - strelka - OK
#CALLER="--caller strelka"
#$MAIN_FILTER --local-script $VAF_FILTER_LOCAL $STRELKA_VCF $SNV_VAF_ARGS $CALLER

