from __future__ import print_function
import sys
import vcf.filters

# Portable printing to stderr, from https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python-2
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class DepthFilter(vcf.filters.Base):
    'Filter variant sites by read depth'
    # Normally we would be able to use the built-in filter "dps"; however, pindel does not write the DP tag and depth filtering fails

    name = 'read_depth'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min_depth', type=int, default=0, help='Retain sites where read depth for tumor and normal > given value')
#        parser.add_argument('--depth_debug', action="store_true", default=False, help='Print debugging information to stderr')
#        parser.add_argument('--depth_caller', type=str, required=True, choices=['strelka', 'varscan', 'pindel'], help='Caller type')

        parser.add_argument('--tumor_name', type=str, default="TUMOR", help='Tumor sample name in VCF')
        parser.add_argument('--normal_name', type=str, default="NORMAL", help='Normal sample name in VCF')
        parser.add_argument('--caller', type=str, required=True, choices=['strelka', 'varscan', 'pindel'], help='Caller type')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')

    def __init__(self, args):
        self.min_depth = args.min_depth
        self.tumor_name = args.tumor_name
        self.normal_name = args.normal_name
        self.caller = args.caller
        self.debug = args.debug

        # below becomes Description field in VCF
        self.__doc__ = "Retain calls where read depth in tumor and normal > %d " % (self.min_depth)

    def filter_name(self):
        return self.name

    def get_depth_strelka(self, VCF_data):
        depth = VCF_data.DP
        if self.debug:
            eprint("strelka depth = %d" % depth)
        return depth

    def get_depth_varscan(self, VCF_data):
        depth = VCF_data.DP
        if self.debug:
            eprint("varscan depth = %d" % depth)
        return depth

    def get_depth_pindel(self, VCF_data):
        rc_ref, rc_var = VCF_data.AD
        depth = rc_ref + rc_var
        if self.debug:
            eprint("pindel depth = %d" % depth)
        return depth

    def get_depth(self, VCF_record, sample_name):
        data=VCF_record.genotype(sample_name).data
        variant_caller = self.caller  
        if self.debug:
            eprint(variant_caller + sample_name)
        if variant_caller == 'strelka':
            return self.get_depth_strelka(data)
        elif variant_caller == 'varscan':
            return self.get_depth_varscan(data)
        elif variant_caller == 'pindel':
            return self.get_depth_pindel(data)
        else:
            raise Exception( "Unknown caller: " + variant_caller)

    def __call__(self, record):

        depth_N = self.get_depth(record, self.normal_name)
        depth_T = self.get_depth(record, self.tumor_name)

        if (self.debug):
            eprint("Variant, Reference depths: %d, %d" % (depth_N, depth_T))

        if depth_N < self.min_depth:
            if (self.debug): eprint("Failed NORMAL min_depth = %d" % depth_N)
            return "depth_N: %d" % depth_N
        if depth_T < self.min_depth:
            if (self.debug): eprint("Failed TUMOR min_depth = %d" % depth_T)
            return "depth_T: %d" % depth_T

        if (self.debug):
            eprint("Passes read depth filter")

