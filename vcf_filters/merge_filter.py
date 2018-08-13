from __future__ import print_function
import sys
import vcf.filters

# filter to include or exclude calls based on their caller, as defined by INFO field "set"

# Portable printing to stderr, from https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python-2
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class MergedCallerFilter(vcf.filters.Base):
    'Filter variant sites by caller, as defined by INFO field "set".'

    name = 'merged_caller'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--include', help='Retain only calls with given caller(s); comma-separated list')
        parser.add_argument('--exclude', help='Exclude all calls with given caller(s); comma-separated list')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')

    def __init__(self, args):

        # User defines either include caller or exclude caller, but not both
        if bool(args.include) == bool(args.exclude):
            raise Exception("Must define exactly one of the following: --include, --exclude")

        if args.include is not None:
            self.including = True
            self.callers = args.include.split(',') 
        else:
            self.including = False
            self.callers = args.exclude.split(',') 

        # below becomes Description field in VCF
        if self.including:
            self.__doc__ = "Retain calls where 'set' INFO field includes one of " + args.include
        else:
            self.__doc__ = "Exclude calls where 'set' INFO field includes any of " + args.exclude

        self.debug = args.debug

    def filter_name(self):
        return self.name

    def __call__(self, record):
        # "caller" is defined by "set" info field
        caller = record.INFO['set']

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
