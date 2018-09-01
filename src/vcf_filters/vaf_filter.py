from common_filter import *
import sys

# Filter VCF files according to tumor, normal VAF values
#
# This is called from pyvcf's `vcf_filter.py` with `vaf` module.
# the following parameters are required:
# * min_vaf_somatic
# * max_vaf_germline
# * tumor_name
# * normal_name
#
# These may be specified on the command line (e.g., --min_vaf_somatic 0.05) or in
# configuration file, as specified by --config config.ini  Sample contents of config file:
#   [vaf]
#   min_vaf_somatic = 0.05
#
# Required command line parameter:
# --caller caller - specifies tool used for variant call. 'strelka', 'varscan', 'pindel', 'merged'
#
# optional command line parameters
# --debug
# --config config.ini
# --bypass

class TumorNormal_VAF(ConfigFileFilter):
    'Filter variant sites by tumor and normal VAF (variant allele frequency)'

    name = 'vaf'

    # for merged caller, will evaluate call based on value of 'set' variable

    @classmethod
    def customize_parser(self, parser):
        # super(TumorNormal_VAF, cls).customize_parser(parser)

        parser.add_argument('--min_vaf_somatic', type=float, help='Retain sites where tumor VAF > than given value')
        parser.add_argument('--max_vaf_germline', type=float, help='Retain sites where normal VAF <= than given value')
        parser.add_argument('--tumor_name', type=str, help='Tumor sample name in VCF')
        parser.add_argument('--normal_name', type=str, help='Normal sample name in VCF')
        parser.add_argument('--caller', type=str, required=True, choices=['strelka', 'varscan', 'pindel', 'merged'], help='Caller type')
        parser.add_argument('--config', type=str, help='Optional configuration file')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')
        parser.add_argument('--bypass', action="store_true", default=False, help='Bypass filter by retaining all variants')
        
    def __init__(self, args):
        # These will not be set from config file (though could be)
        self.caller = args.caller
        self.debug = args.debug
        self.bypass = args.bypass

        # Read arguments from config file first, if present.
        # Then read from command line args, if defined
        # Note that default values in command line args would
        #   clobber configuration file values so are not defined
        config = self.read_config_file(args.config)

        self.set_args(config, args, "min_vaf_somatic", arg_type="float")
        self.set_args(config, args, "max_vaf_germline", arg_type="float")
        self.set_args(config, args, "tumor_name")
        self.set_args(config, args, "normal_name")

        # below becomes Description field in VCF
        if self.bypass:
            self.__doc__ = "Bypassing Tumor Normal VAF filter, retaining all reads"
        else:
            self.__doc__ = "Retain calls where normal VAF <= %f and tumor VAF >= %f " % (self.max_vaf_germline, self.min_vaf_somatic)
            
    def filter_name(self):
        return self.name

    def get_readcounts_strelka(self, VCF_record, VCF_data):
        # pass VCF_record only to extract info (like ALT and is_snp) not available in VCF_data

        if not VCF_record.is_snp:
            raise Exception( "Only SNP calls supported for Strelka: " + VCF_record)
        # read counts, as dictionary. e.g. {'A': 0, 'C': 0, 'T': 0, 'G': 25}
        tier=0   
        rc = {'A':VCF_data.AU[tier], 'C':VCF_data.CU[tier], 'G':VCF_data.GU[tier], 'T':VCF_data.TU[tier]}

        # Sum read counts across all variants. In some cases, multiple variants are separated by , in ALT field
        # Implicitly, only SNV supported here.
        #   Note we convert vcf.model._Substitution to its string representation to use as key
        rc_var = sum( [rc[v] for v in map(str, VCF_record.ALT) ] )
        rc_tot = sum(rc.values())
        vaf = float(rc_var) / float(rc_tot)
        if self.debug:
            eprint("rc: %s, rc_var: %f, rc_tot: %f, vaf: %f" % (str(rc), rc_var, rc_tot, vaf))
        return vaf

    def get_readcounts_varscan(self, VCF_record, VCF_data):
        # We'll take advantage of pre-calculated VAF
        # Varscan: CallData(GT=0/0, GQ=None, DP=96, RD=92, AD=1, FREQ=1.08%, DP4=68,24,1,0)
        ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
        # This works for both snp and indel calls
        vaf = VCF_data.FREQ
        if self.debug:
            eprint("VCF_data.FREQ = %s" % vaf)
        return float(vaf.strip('%'))/100.

    def get_readcounts_pindel(self, VCF_record, VCF_data):
        # read counts supporting reference, variant, resp.
        rc_ref, rc_var = VCF_data.AD
        vaf = rc_var / float(rc_var + rc_ref)
        if self.debug:
            eprint("pindel VCF = %f" % vaf)
        return vaf

    def get_vaf(self, VCF_record, sample_name, variant_caller=None):
        data=VCF_record.genotype(sample_name).data
        if variant_caller is None:
            variant_caller = self.caller  # we permit the possibility that each line has a different caller

        if variant_caller == 'strelka':
            return self.get_readcounts_strelka(VCF_record, data)
        elif variant_caller == 'varscan':
            return self.get_readcounts_varscan(VCF_record, data)
        elif variant_caller == 'pindel':
            return self.get_readcounts_pindel(VCF_record, data)
        elif variant_caller == 'merged':
            # Caller is contained in 'set' INFO field
            merged_caller = VCF_record.INFO['set'][0]
            # TODO: finish and test this thoroughly.  In a merged line with set=A-B-C, call readcounts with "A"
            if merged_caller == 'strelka':
                return self.get_readcounts_strelka(VCF_record, data)
            elif merged_caller == 'varscan':
                return self.get_readcounts_varscan(VCF_record, data)
            elif merged_caller == 'varindel':
                return self.get_readcounts_varscan(VCF_record, data)
            elif merged_caller == 'strelka-varscan':
                return self.get_readcounts_strelka(VCF_record, data)
            elif merged_caller == 'pindel':
                return self.get_readcounts_pindel(VCF_record, data)
            else:
                raise Exception( "Unknown caller in INFO set field: " + merged_caller)
        else:
            raise Exception( "Unknown caller: " + variant_caller)

    def __call__(self, record):
        vaf_N = self.get_vaf(record, self.normal_name)
        vaf_T = self.get_vaf(record, self.tumor_name)

        if (self.debug):
            eprint("Normal, Tumor vaf: %f, %f" % (vaf_N, vaf_T))

        if self.bypass:
            if (self.debug): eprint("** Bypassing filter, retaining read **" )
            return

##       Original logic, with 2=Tumor
##           RETAIN if($rc2var/$r_tot2>=$min_vaf_somatic && $rcvar/$r_tot<=$max_vaf_germline && $r_tot2>=$min_coverage && $r_tot>=$min_coverage)
##       Here, logic is reversed.  We return if fail a test
        if vaf_T < self.min_vaf_somatic:
            if (self.debug):
                eprint("** Failed vaf_T < min_vaf_somatic **")
            return "VAF_T: %f" % vaf_T
        if vaf_N >= self.max_vaf_germline:
            if (self.debug):
                eprint("** Failed vaf_N >= max_vaf_germline **")
            return "VAF_N: %f" % vaf_N
        if (self.debug):
            eprint("** Passes VAF filter **")

