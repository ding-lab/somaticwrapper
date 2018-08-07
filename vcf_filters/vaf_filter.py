from __future__ import print_function
import sys
import vcf.filters

# Portable printing to stderr, from https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python-2
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class TumorNormal_VAF(vcf.filters.Base):
    'Filter variant sites by tumor and normal VAF (variant allele frequency)'

    name = 'vaf'

##       RETAIN if($rcTvar/$r_totT>=$min_vaf_somatic && $rcvar/$r_tot<=$max_vaf_germline && $r_totT>=$min_coverage && $r_tot>=$min_coverage)

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min_vaf_somatic', type=float, default=0.05, help='Retain sites where tumor VAF > than given value')
        parser.add_argument('--max_vaf_germline', type=float, default=0.02, help='Retain sites where normal VAF <= than given value')
        parser.add_argument('--tumor_name', type=str, default="TUMOR", help='Tumor sample name in VCF')
        parser.add_argument('--normal_name', type=str, default="NORMAL", help='Normal sample name in VCF')
        parser.add_argument('--caller', type=str, required=True, choices=['strelka', 'varscan', 'pindel'], help='Caller type')
        parser.add_argument('--debug', action="store_true", default=False, help='Print debugging information to stderr')

    def __init__(self, args):
        self.min_vaf_somatic = args.min_vaf_somatic
        self.max_vaf_germline = args.max_vaf_germline
        self.tumor_name = args.tumor_name
        self.normal_name = args.normal_name
        # below becomes Description field in VCF
        self.__doc__ = "Retain calls where normal VAF <= %f and tumor VAF >= %f " % (self.max_vaf_germline, self.min_vaf_somatic)
        self.caller = args.caller
        self.debug = args.debug
            
    def filter_name(self):
        return self.name

    def get_readcounts_strelka(self, VCF_data):
        # read counts, as dictionary. e.g. {'A': 0, 'C': 0, 'T': 0, 'G': 25}
        tier=0   
        rc = {'A':VCF_data.AU[tier], 'C':VCF_data.CU[tier], 'G':VCF_data.GU[tier], 'T':VCF_data.TU[tier]}

#        TODO: confirm this is a SNV

        # Sum read counts across all variants. In some cases, multiple variants are separated by , in ALT field
        # Implicitly, only SNV supported here.
        #   Note we convert vcf.model._Substitution to its string representation to use as key
        rc_var = sum( [rc[v] for v in map(str, VCF_record.ALT) ] )
        rc_tot = sum(rc.values())
        vaf = float(rc_var) / float(rc_tot)
        if self.debug:
            eprint("rc: %s, rc_var: %f, rc_tot: %f, vaf: %f" % (str(rc), rc_var, rc_tot, vaf))
        return vaf

    def get_readcounts_varscan(self, VCF_data):
        # We'll take advantage of pre-calculated VAF
        # Varscan: CallData(GT=0/0, GQ=None, DP=96, RD=92, AD=1, FREQ=1.08%, DP4=68,24,1,0)
        ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
        # This works for both snp and indel calls
        vaf = VCF_data.FREQ
        if self.debug:
            eprint("VCF_data.FREQ = %s" % vaf)
        return float(vaf.strip('%'))/100.

    def get_readcounts_pindel(self, VCF_data):
        # read counts supporting reference, variant, resp.
        rc_ref, rc_var = VCF_data.AD
        vaf = rc_var / float(rc_var + rc_ref)
        if self.debug:
            eprint("pindel VCF = %f" % vaf)
        return vaf

    def get_vaf(self, VCF_record, sample_name):
        data=VCF_record.genotype(sample_name).data
        variant_caller = self.caller  # we permit the possibility that each line has a different caller
        if self.debug:
            eprint(variant_caller + sample_name)
        if variant_caller == 'strelka':
            if VCF_record.is_snp() is False:
                raise Exception( "Only SNP calls supported for Strelka: " + VCF_record)
            return self.get_readcounts_strelka(data)
        elif variant_caller == 'varscan':
            return self.get_readcounts_varscan(data)
        elif variant_caller == 'pindel':
            return self.get_readcounts_pindel(data)
        else:
            raise Exception( "Unknown caller: " + variant_caller)

    def __call__(self, record):
        vaf_N = self.get_vaf(record, self.normal_name)
        vaf_T = self.get_vaf(record, self.tumor_name)

##       Original logic, with 2=Tumor
##       RETAIN if($rc2var/$r_tot2>=$min_vaf_somatic && $rcvar/$r_tot<=$max_vaf_germline && $r_tot2>=$min_coverage && $r_tot>=$min_coverage)
##       Here, logic is reversed.  We return if fail a test
        if vaf_T < self.min_vaf_somatic:
            if (self.debug):
                eprint("Failed vaf_T < min_vaf_somatic")
            return "VAF_T: %f" % vaf_T
        if vaf_N >= self.max_vaf_germline:
            if (self.debug):
                eprint("Failed vaf_N >= max_vaf_germline")
            return "VAF_N: %f" % vaf_N
        if (self.debug):
            eprint("Passes VAF filter")

class IndelLengthFilter(vcf.filters.Base):
    'Filter indel variant sites by reference and variant length'

    name = 'indel_length'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min_length', type=int, default=0, help='Retain sites where indel length > given value')
        parser.add_argument('--max_length', type=int, default=0, help='Retain sites where indel length <= given value. 0 disables test')
        parser.add_argument('--length_debug', action="store_true", default=False, help='Print debugging information to stderr')

    def __init__(self, args):
        self.min_length = args.min_length
        self.max_length = args.max_length
        # below becomes Description field in VCF
        self.__doc__ = "Retain calls where indel length > %d and < %d " % (self.min_length, self.max_length)
        self.debug = args.length_debug

    def filter_name(self):
        return self.name

    def __call__(self, record):

        len_ALT = len(record.ALT)
        len_REF = len(record.REF)
        if (self.debug):
            eprint("Variant, Reference: %s, %s" % (record.ALT, record.REF))
            eprint("Variant, Reference lengths: %d, %d" % (len_ALT, len_REF))

        if len_ALT < self.min_length:
            if (self.debug): eprint("Failed ALT min_length = %d" % len_ALT)
            return "len_ALT: %f" % len_ALT
        if len_REF < self.min_length:
            if (self.debug): eprint("Failed REF min_length = %d" % len_REF)
            return "len_REF: %f" % len_REF

        if self.max_length is not 0:
            if len_ALT > self.max_length:
                if (self.debug): eprint("Failed ALT max_length = %d" % len_ALT)
                return "len_ALT: %f" % len_ALT
            if len_REF > self.max_length:
                if (self.debug): eprint("Failed REF max_length = %d" % len_REF)
                return "len_REF: %f" % len_REF

        if (self.debug):
            eprint("Passes length filter")

class DepthFilter(vcf.filters.Base):
    'Filter variant sites by read depth'
    # Normally we would be able to use the built-in filter "dps"; however, pindel does not write the DP tag and depth filtering fails

    name = 'read_depth'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min_depth', type=int, default=0, help='Retain sites where read depth for tumor and normal > given value')
        parser.add_argument('--depth_debug', action="store_true", default=False, help='Print debugging information to stderr')
        parser.add_argument('--depth_caller', type=str, required=True, choices=['strelka', 'varscan', 'pindel'], help='Caller type')

    def __init__(self, args):
        self.min_depth = args.min_depth
        # below becomes Description field in VCF
        self.__doc__ = "Retain calls where read depth in tumor and normal > %d " % (self.min_depth)
        self.debug = args.depth_debug

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

