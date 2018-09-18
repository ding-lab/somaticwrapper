from common_filter import *
import sys

# filter to include or exclude calls based on their caller, as defined by INFO field "set"

class MergeFilter(ConfigFileFilter):
    'Filter variant sites by caller, as defined by INFO field "set".'

    name = 'merge'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--include', help='Retain only calls with given caller(s); comma-separated list')
        parser.add_argument('--exclude', help='Exclude all calls with given caller(s); comma-separated list')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')
        parser.add_argument('--bypass', action="store_true", default=False, help='Bypass filter by retaining all variants')

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

        self.debug = args.debug
        self.bypass = args.bypass

        # below becomes Description field in VCF
        if self.bypass:
            self.__doc__ = "Bypassing Merge filter, retaining all reads"
        elif self.including:
            self.__doc__ = "Retain calls where 'set' INFO field includes one of " + args.include
        else:
            self.__doc__ = "Exclude calls where 'set' INFO field includes any of " + args.exclude

    def filter_name(self):
        return self.name

    def __call__(self, record):
        # "caller" is defined by "set" info field
        caller = record.INFO['set']

        if self.bypass:
            if (self.debug): eprint("** Bypassing %s filter, retaining read **" % self.name )
            return

        if self.including:
            # keep call only if caller is in callers list
            if caller not in self.callers:
                if (self.debug): eprint("** FAIL: %s not in %s **" % (caller, str(self.callers)))
                return "unknown " + caller
            else:
                if (self.debug): eprint("** PASS: %s in %s **" % (caller, str(self.callers)))
                return
        else:                 
            # keep call only if caller is not in callers list
            if caller in self.callers:
                if (self.debug): eprint("** FAIL: %s is in %s **" % (caller, str(self.callers)))
                return "excluding " + caller
            else:
                if (self.debug): eprint("** PASS: %s not in %s **" % (caller, str(self.callers)))
                return
