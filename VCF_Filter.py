#!/usr/bin/python

# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu

# Filter VCF files by variant allele frequency (VAF) and length


import vcf # http://pyvcf.readthedocs.io/en/latest/API.html
import sys

def write_bpc(o, vcf_reader, attribute, no_header=False):
    """Write BreakpointSurveyor BPC file"""
        
    for record in vcf_reader:
        chromA, posA, chromB, posB = record.CHROM, record.start, record.INFO['CHR2'], record.sv_end


        o.write('\t'.join( map(str, linedata) ) + "\n")

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] mode ...
        """

    parser = OptionParser(usage_text, version="$Revision: 1.2 $")
    parser.add_option("-i", dest="infn", default="stdin", help="Input filename")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output filename")

    (options, params) = parser.parse_args()

#    if (len(params) != 1):
#        parser.error("Pass 1 argument.")
#    mode = params[0]
#    if mode not in ['bed', 'size', 'bpc']:
#        parser.error("mode must be one of ('bed', 'size', 'bpc')")

    if options.infn == "stdin":
        f = sys.stdin
    else:
        f = open(options.infn, 'r')
    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    vcf_reader = vcf.Reader(f)

    filter_vcf(f, o)

    f.close()
    o.close()

if __name__ == '__main__':
    main()

