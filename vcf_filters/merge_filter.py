from __future__ import print_function
import sys
import vcf.filters

# filter to include or exclude calls based on their caller, as defined by INFO field "set"

# Portable printing to stderr, from https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python-2
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class IndelLengthFilter(vcf.filters.Base):
    'Filter variant sites by caller, as defined by INFO field "set".'

    name = 'indel_length'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--include_caller', help='Retain only calls with given caller(s); comma-separated list')
        parser.add_argument('--exclude_caller', help='Exclude all calls with given caller(s); comma-separated list')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')

    def __init__(self, args):
        # define either include_caller or exclude_caller; test with XOR        
        if (args.include_caller is None) ^ (args.exclude_caller is None):
            raise Exception("Must define exactly one of the following: --include_caller, --exclude_caller")

        if args.include_caller is not None:
            self.including = True
            self.callers = self.include_caller.split(',') 
        else:
            self.including = False
            self.callers = self.exclude_caller.split(',') 

        # below becomes Description field in VCF
        if self.including:
            self.__doc__ = "Retain calls where 'set' INFO field includes one of " + args.include_caller
        else:
            self.__doc__ = "Exclude calls where 'set' INFO field includes any of " + args.exclude_caller

        self.debug = args.debug

    def filter_name(self):
        return self.name

    def __call__(self, record):
        # "caller" is defined by "set" info field
        assert len(VCF_record.INFO['set']) == 1 # assuming that has only one field, need to rework comparison logic if this breaks
        caller = VCF_record.INFO['set'][0]

        if self.including:
            # keep call only if caller is in callers list
            if caller not in self.callers:
                if (self.debug): eprint("** Failed: %s not in %s **" % (caller, str(self.callers)))
                return "unknown " + caller
            else:
                if (self.debug): eprint("** Passes: %s in %s **" % (caller, str(self.callers)))
                return
        else:                 
            # keep call only if caller is not in callers list
            if caller in self.callers:
                if (self.debug): eprint("** Failed: %s is in %s **" % (caller, str(self.callers)))
                return "excluding " + caller
            else:
                if (self.debug): eprint("** Passes: %s not in %s **" % (caller, str(self.callers)))
                return
