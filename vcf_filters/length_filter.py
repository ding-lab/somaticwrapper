from common_filter import *
import sys

# Filter VCF files according to indel length:
# Retain calls where length of variant and reference > min_length and < max_length 
#
# The following parameters are required:
# * min_length
# * max_length
#
# These may be specified on the command line (e.g., --max_length 100) or in
# configuration file, as specified by --config config.ini  Sample contents of config file:
#   [indel_length]
#   max_length = 100
#
# optional command line parameters
# --debug
# --config config.ini

class IndelLengthFilter(ConfigFileFilter):
    'Filter indel variant sites by reference and variant length'

    name = 'indel_length'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min_length', type=int, help='Retain sites where indel length > given value')
        parser.add_argument('--max_length', type=int, help='Retain sites where indel length <= given value. 0 disables test')
        parser.add_argument('--config', type=str, help='Optional configuration file')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')

    def __init__(self, args):
        # These will not be set from config file (though could be)
        self.debug = args.debug

        # Read arguments from config file first, if present.
        # Then read from command line args, if defined
        # Note that default values in command line args would
        #   clobber configuration file values so are not defined

        config = self.read_config_file(args.config)

        self.set_args(config, args, "min_length", arg_type="int")
        self.set_args(config, args, "max_length", arg_type="int")

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
            if (self.debug): eprint("** Failed REF min_length = %d **" % len_REF)
            return "len_REF: %f" % len_REF
        if len_ALT < self.min_length:
            if (self.debug): eprint("** Failed ALT min_length = %d **" % len_ALT)
            return "len_ALT: %f" % len_ALT

        if self.max_length is not 0:
            if len_REF > self.max_length:
                if (self.debug): eprint("** Failed REF max_length = %d **" % len_REF)
                return "len_REF: %f" % len_REF
            if len_ALT > self.max_length:
                if (self.debug): eprint("** Failed ALT max_length = %d **" % len_ALT)
                return "len_ALT: %f" % len_ALT

        if (self.debug):
            eprint("** Passes length filter **")

