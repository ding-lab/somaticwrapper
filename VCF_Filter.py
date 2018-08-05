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


# Pindel after merge
# 1   14106395    .   C   CTCA    .   PASS    AC=0;AF=0.00;AN=4;END=14106395;HOMLEN=2;HOMSEQ=TC;SVLEN=3;SVTYPE=INS;set=pindel GT:AD   ./. ./. 0/0:11,0    0/0:16,2
# Varscan after merge 
# 1   16911827    .   C   G   .   PASS    AC=1;AF=0.250;AN=4;DP=267;GPV=1E0;SOMATIC;SPV=2.0555E-2;SS=2;SSC=16;set=varscan GT:AD:DP:DP4:FREQ:RD    0/0:2:109:105,2,2,0:1.83%:107   0/1:13:158:140,5,13,0:8.23%:145 ./. ./.


def filter_vcf(o, vcf_reader, no_header=False):
    '''Write BreakpointSurveyor BPC file'''
        
    for record in vcf_reader:

        # rc_N, rc_T are tuples with read counts for tier 1 A, C, G, T for Normal and Tumor, resp
        data_N=record.genotype('NORMAL').data
        data_T=record.genotype('TUMOR').data
        rc_N = {'A':data_N.AU[0], 'C':data_N.CU[0], 'G':data_N.GU[0], 'T':data_N.TU[0]}
        rc_T = {'A':data_T.AU[0], 'C':data_T.CU[0], 'G':data_T.GU[0], 'T':data_T.TU[0]}

        # Sum all variants. Need to convert vcf.model._Substitution to its string representation
        rc_var_N = sum( [rc_N[v] for v in map(str, record.ALT) ] )
        rc_var_T = sum( [rc_T[v] for v in map(str, record.ALT) ] )

        eprint(rc_var_N)
        eprint(rc_var_T)
        sys.exit(1)


        for v in map(str, record.ALT): # Convert vcf.model._Substitution to its string representation
            if v not in rc_N:
                eprint("Warning: rc_N does not have counts for variant %s" % v)
            eprint(rc_N[v])
            eprint(rc_T[v])
#            eprint(rc_N[v])
#            if v not in rc_T.keys():
#                eprint("Warning: rc_T does not have counts for variant %s" % v)
#            eprint(rc_T[v])

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

