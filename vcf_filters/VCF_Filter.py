#!/usr/bin/python

# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu

# Reproduce functionality of SomaticWrapper:master script vaf_filter_v1.2.pl
# https://github.com/ding-lab/somaticwrapper/blob/master/vaf_filter_v1.2.pl

# Filter VCF files by 
# * variant allele frequency (VAF) 
# * indel length
# * Coverage


from __future__ import print_function
import vcf # http://pyvcf.readthedocs.io/en/latest/API.html
import sys
import argparse # https://docs.python.org/2/library/argparse.html

# Portable printing to stderr, from https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python-2
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# Returns tuple of normal, tumor VCF
# VAF: (number of reads supporting variant)/(number of reads supporting variant and number of reads supporting reference)
def get_vaf(record):
    #chromA, posA, chromB, posB = record.CHROM, record.start, record.INFO['CHR2'], record.sv_end
    pass


def get_readcounts_SNV(VCF_call):
```
For given call, return (read_counts, read_counts_variant, read_counts_total),
where read_counts is dictionary of read counts by nucleotide, while read_counts_variant and read_counts_total
are sums of counts of variant (from record.ALT) 
pyvcf Call: http://pyvcf.readthedocs.io/en/latest/API.html#vcf-model-call
    Can be obtained with record.genotype(name)
```
    # rc is dictionary with read counts for tier1 A, C, G, T 
    tier=0   # this may be specific to STRELKA.  Want tier1
    data=VCF_call.data

    # read counts, as dictionary. e.g. {'A': 0, 'C': 0, 'T': 0, 'G': 25}
    rc = {'A':data.AU[tier], 'C':data.CU[tier], 'G':data.GU[tier], 'T':data.TU[tier]}

    # Sum read counts across all variants. In some cases, multiple variants are separated by , in ALT field
    # Implicitly, only SNV supported here.
    #   Note we convert vcf.model._Substitution to its string representation to use as key
    rc_var = sum( [rc[v] for v in map(str, record.ALT) ] )
    rc_tot = sum(d.values(rc))

    return (rc, rc_var, rc_tot)


def write_stats
print OUT2 $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t",$info,"\t",$rc{$ref},"\t",$rc{$ref}/$r_tot,"\t",$rcvar,"\t",$rcvar/$r_tot,"\t",$rc2{$ref},"\t",$rc2{$ref}/$r_tot2,"\t",$rc2var,"\t",$rc2var/$r_tot2,"\n";



# Modeled after ~/strelka-varscan/ part of vaf_filter_v1.2.pl
# Pindel after merge
# 1   14106395    .   C   CTCA    .   PASS    AC=0;AF=0.00;AN=4;END=14106395;HOMLEN=2;HOMSEQ=TC;SVLEN=3;SVTYPE=INS;set=pindel GT:AD   ./. ./. 0/0:11,0    0/0:16,2
# Varscan after merge 
# 1   16911827    .   C   G   .   PASS    AC=1;AF=0.250;AN=4;DP=267;GPV=1E0;SOMATIC;SPV=2.0555E-2;SS=2;SSC=16;set=varscan GT:AD:DP:DP4:FREQ:RD    0/0:2:109:105,2,2,0:1.83%:107   0/1:13:158:140,5,13,0:8.23%:145 ./. ./.


def filter_vcf(o, vcf_reader, no_header=False):
    '''Write BreakpointSurveyor BPC file'''

        
    for record in vcf_reader:
        VCF_header = record.CHROM, record.POS, record.REF, record.ALT, record.FILTER

        # Get read counts for tumor and normal samples
        rc_N, rc_var_N, rc_tot_N = get_readcounts_snv(record.genotype('NORMAL'))
        rc_T, rc_var_T, rc_tot_T = get_readcounts_snv(record.genotype('TUMOR'))

### calculate vaf for tumor and normal:
## vaf calculation: defined as (number of reads supporting variant)/(number of reads supporting variant and number of reads supporting reference) ####
## use 20 a coverage filtering ##
        vaf_N = rc_var_N / rc_tot_N
        vaf_T = rc_var_T / rc_tot_T

#            if($rc2var/$r_tot2>=$min_vaf_somatic && $rcvar/$r_tot<=$max_vaf_germline && $r_tot2>=$min_coverage && $r_tot>=$min_coverage) 
#            {
#                print OUT1 $ltr,"\n";   # merged.filtered.vcf 
#            } 	

        sys.exit(1)

        #o.write('\t'.join( map(str, linedata) ) + "\n")
        eprint("Hello")
        eprint(record.INFO)
#        eprint(record.FORMAT)
        eprint(record.samples)
        eprint(record.genotype('NORMAL').data.GU)
        eprint(record.genotype('TUMOR').data.GU)
#        eprint(record.samples[0].data.DP)
        sys.exit(1)
#        eprint("Normal VAF: %f  Tumor VAF: %f" % get_vcf(record)

def main():
    usage_text = '''
    Test parser for VCFs
    '''

    parser = argparse.ArgumentParser( description=usage_text ) # See also https://pymotw.com/2/argparse/
    parser.add_argument( '--version', action='version', version='%(prog)s 1.0' )
    parser.add_argument( '-i', dest='infn', default='stdout', help='Input VCF (default reads stdin)' )
    parser.add_argument( '-o', dest='outfn', default='stdout', help='Output Pass VCF (default writes to stdout)' )
    parser.add_argument( '-H', action='store_true', dest='noWriteHeader', default=False, help='Do not write header to output' )
#    parser.add_argument( dest='ROI', metavar='ROI_file', type=str, help='population covered ROI file name' )
#    parser.add_argument( dest='vcf', metavar='VCF_file', type=str, help='VCF file containing mutation calls for all the samples in the population' )

    args = parser.parse_args()

#    if (len(params) != 1):
#        parser.error("Pass 1 argument.")
#    mode = params[0]
#    if mode not in ['bed', 'size', 'bpc']:
#        parser.error("mode must be one of ('bed', 'size', 'bpc')")

    if args.infn == "stdin":
        f = sys.stdin
    else:
        f = open(args.infn, 'r')
    if args.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(args.outfn, "w")

    vcf_reader = vcf.Reader(f)

    filter_vcf(o, vcf_reader, args.noWriteHeader)

    f.close()
    o.close()

if __name__ == '__main__':
    main()

