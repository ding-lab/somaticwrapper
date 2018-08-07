from __future__ import print_function
import sys
import vcf.filters

# Portable printing to stderr, from https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python-2
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class IndelLengthFilter(vcf.filters.Base):
    'Filter indel variant sites by reference and variant length'

    name = 'indel_length'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min_length', type=int, default=0, help='Retain sites where indel length > given value')
        parser.add_argument('--max_length', type=int, default=0, help='Retain sites where indel length <= given value. 0 disables test')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')

    def __init__(self, args):
        self.min_length = args.min_length
        self.max_length = args.max_length
        # below becomes Description field in VCF
        self.__doc__ = "Retain calls where indel length > %d and < %d " % (self.min_length, self.max_length)
        self.debug = args.debug

    def filter_name(self):
        return self.name

    def __call__(self, record):

        # ALT will typically be a list
        if isinstance(record.ALT, list):
            len_ALT = len(max(record.ALT, key=len)) 
        else:
            len_ALT = len(record.ALT)
        len_REF = len(record.REF)
        if (self.debug):
            eprint("Reference, Variant: %s, %s" % (record.REF, record.ALT))
            eprint("Reference, Variant lengths: %d, %d" % (len_REF, len_ALT))

        if len_REF < self.min_length:
            if (self.debug): eprint("Failed REF min_length = %d" % len_REF)
            return "len_REF: %f" % len_REF
        if len_ALT < self.min_length:
            if (self.debug): eprint("Failed ALT min_length = %d" % len_ALT)
            return "len_ALT: %f" % len_ALT

        if self.max_length is not 0:
            if len_REF > self.max_length:
                if (self.debug): eprint("Failed REF max_length = %d" % len_REF)
                return "len_REF: %f" % len_REF
            if len_ALT > self.max_length:
                if (self.debug): eprint("Failed ALT max_length = %d" % len_ALT)
                return "len_ALT: %f" % len_ALT

        if (self.debug):
            eprint("Passes length filter")

