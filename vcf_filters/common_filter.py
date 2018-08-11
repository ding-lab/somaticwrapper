from __future__ import print_function
import sys
import vcf.filters
import ConfigParser
import os.path


def eprint(*args, **kwargs):
# Portable printing to stderr, from https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python-2
    print(*args, file=sys.stderr, **kwargs)

# Configuration file and initialization of parameters:
# * we can read a configuration file (ini format) as an alternative to passing command line parameters
#   * --config config.ini will read config file config.ini, and read parameters associated with section [ vaf ]
#   * See details here: https://docs.python.org/3/library/configparser.html#supported-ini-file-structure
# * A parameter in configuration file will be overridden by command line argument of the same name
# * No default command line values in general

class ConfigFileFilter(vcf.filters.Base):
    'Base class of pyvcf filters which can read from configuration ini file'

    def set_args(self, config, args, option, required=True, arg_type="string"):
        '''
        Set class attributes directly from command line args (priority) or configparser, if present.
        set_args(self,config,args,"foo") will define self.foo = "foo value"
        arg_type = "float" will cast to float.  This should be generalized.
        '''
        value = None
        if option in vars(args) and vars(args)[option] is not None:
            value = vars(args)[option]
            if self.debug:
                eprint("Setting %s = %s from args" % (option, value))
        elif config is not None and config.has_option(self.name, option):
            if arg_type == "float":
                value = config.getfloat(self.name, option)
            else:
                value = config.get(self.name, option)
            if self.debug:
                eprint("Setting %s = %s from config" % (option, value))

        if value is None:
            msg = "Argument %s not defined" % option
            if required:
                raise Exception("Error: %s " % msg)
            else:
                eprint("Config value %s not defined" % option)
        else:
            setattr(self, option, value)

# https://docs.python.org/3/library/configparser.html and https://docs.python.org/2/library/configparser.html
    def read_config_file(self, config_fn):
    # return None if not defined

        if config_fn is None:
            return None
        eprint("Reading configuration file " + config_fn)
        if not os.path.isfile(config_fn):
            raise Exception("Error: Configuration file %s not found." % config_fn)
        config = ConfigParser.ConfigParser()
        config.read(config_fn)
        return config
