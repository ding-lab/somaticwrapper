import csv
from dataclasses import dataclass
import logging
from pathlib import Path
from textwrap import dedent
from typing import Dict

_logger = logging.getLogger(__name__)

# List of sample names to run the pipeline
SAMPLE_LIST_PTH = config['sample_list']
# The mapping of sample name to local file locations
FILE_MAP_PTH = config['file_map']
WORKFLOW_ROOT = config['workflow_root']  # Path to this repository
BWA_INDEX_PREFIX = config['bwa_index_prefix']  # Path to the BWA index prefix

wildcard_constraints:
    sample="[^/]+"


@dataclass
class SampleInfo:
    """Class to keep track of the sample info."""
    sample_name: str
    case_id: str
    disease: str
    experimental_strategy: str  # WXS or WGS
    gdc_catalog_sample_type: str    # T, N, A, or T.RRNA_rk3EEVp
    sample_type: str    # tumor, blood_normal, tissue_normal
    raw_bam_uuid: str
    raw_bam_pth: Path


# Construct the sample list
SAMPLES = []
with open(SAMPLE_LIST_PTH) as f:
    SAMPLES = [line.strip() for line in f]


# Read file map
SAMPLE_INFO: Dict[str, SampleInfo] = {}
with open(FILE_MAP_PTH) as f:
   reader = csv.DictReader(f, dialect='excel-tab')
   for row in reader:
        sample_name = row['# sample_name']
        case_id = row['case']
        disease = row['disease']
        experimental_strategy = row['experimental_strategy']
        sample_type = row['sample_type']

        # Parse Matt's sample name
        if m := re.search(r"W[XG]S\.([\.\w]+)$", sample_name):
            gdc_catalog_sample_type = m.group(1)
        else:
            raise ValueError(f"Could not parse sample name for {sample_name}")

        raw_bam_uuid = row['UUID']
        raw_bam_pth = Path(row['data_path'])
	print(raw_bam_pth)
	assert raw_bam_pth.exists()

        SAMPLE_INFO[sample_name] = SampleInfo(
            sample_name=sample_name,
            case_id=case_id,
            disease=disease,
            experimental_strategy=experimental_strategy,
            gdc_catalog_sample_type=gdc_catalog_sample_type,
            sample_type=sample_type,
            raw_bam_uuid=raw_bam_uuid,
            raw_bam_pth=raw_bam_pth,
        )


def find_sample_bam(wildcards):
    """Find the BAM file path of a given sample."""
    sample_info = SAMPLE_INFO[wildcards.sample]
    return {
        'bam': str(sample_info.raw_bam_pth),
    }


checkpoint bam_to_fastqs:
    """Export reads as FASTQs from an aligned BAM."""
    output: sample_folder=directory('bam_to_fastqs/{sample}')
    input: unpack(find_sample_bam)
    log: 'logs/bam_to_fastqs/{sample}.log'
    resources:
        io_heavy=1
    shell:
        'mkdir {output.sample_folder}; '
        'bamtofastq '
        'filename={input.bam} inputformat=bam '
        'collate=1 '
        'exclude=QCFAIL,SECONDARY,SUPPLEMENTARY '
        'gz=1 level=5 '
        'tryoq=1 '
        'combs=1 '
        'T=$(mktemp --tmpdir tmp_bamtofastq.XXXXXXXXXXXXXXXX) '  # Write temporary files under $TMPDIR
        'outputdir={output.sample_folder} '
        'outputperreadgroup=1 '
        'outputperreadgroupsuffixF=_1.fq.gz '
        'outputperreadgroupsuffixF2=_2.fq.gz '
        'outputperreadgroupsuffixO=_o1.fq.gz '
        'outputperreadgroupsuffixO2=_o2.fq.gz '
        'outputperreadgroupsuffixS=_s.fq.gz '
        '2>{log} 1>&2'


def calc_required_memory(wildcards, attempt):
    return 16000 + 8000 * (attempt - 1)


rule bwa_align_one_readgroup_fastq:
    """BWA MEM alignemnt of one readgroup FASTQs."""
    output: bam=temp('readgroup_bam/{sample}/{rg}.bam')
    input:
        r1_fq='bam_to_fastqs/{sample}/{rg}_1.fq.gz',
        r2_fq='bam_to_fastqs/{sample}/{rg}_2.fq.gz'
    threads: 8
    resources:
        mem_mb=calc_required_memory
    log: 'logs/bwa_align/{sample}/{rg}.log'
    shell:
        dedent("""
        bwa mem \
            -t {threads} -T 0 \
            -R '@RG\\tID:{wildcards.rg}\\tSM:{wildcards.sample}' \
            {BWA_INDEX_PREFIX} {input.r1_fq} {input.r2_fq} \
            2>{log} \
        | samtools view -Shb -o {output.bam} -
        """)


rule picard_sort_one_readgroup_bam:
    """Picard sort one readgroup BAM."""
    output: bam=temp('readgroup_bam/{sample}/{rg}.sorted.bam')
    input: bam=rules.bwa_align_one_readgroup_fastq.output['bam']
    threads: 4
    resources:
        mem_mb=calc_required_memory
    log: 'logs/picard_sort_rg_bam/{sample}/{rg}.log'
    shell:
        "picard -Xmx{resources.mem_mb}m SortSam "
        "CREATE_INDEX=true "
        "I={input.bam} "
        "O={output.bam} "
        "SORT_ORDER=coordinate "
        "VALIDATION_STRINGENCY=STRICT "
        "TMP_DIR=$(mktemp -d --tmpdir tmp_picard.XXXXXXXXXXXXXXXX) "
        "2>{log} 1>&2"


def get_readgroups_of_one_sample(wildcards):
    """Get the available read groups of a sample."""
    sample_folder = Path(checkpoints.bam_to_fastqs.get(**wildcards).output['sample_folder'])
    all_readgroups = set(p.name[:-len('_1.fq.gz')] for p in sample_folder.glob('*_1.fq.gz'))
    # ignore the `default` RG
    all_readgroups.discard('default')
    return sorted(all_readgroups)


def get_sorted_readgroup_bams_of_one_sample(wildcards):
    all_readgroups = get_readgroups_of_one_sample(wildcards)
    all_sorted_bams = [
        f'readgroup_bam/{wildcards.sample}/{rg}.sorted.bam'
        for rg in all_readgroups
    ]
    return all_sorted_bams


rule picard_merge_readgroup_bams_of_one_sample:
    """Merge all readgroup BAMs of one sample by Picard."""
    output:
        merged_bam=temp('merged_bam/{sample}.bam'),
        merged_bai=temp('merged_bam/{sample}.bai')
    input: sorted_readgroup_bams=get_sorted_readgroup_bams_of_one_sample
    params:
        input_bams_args=lambda wildcards, input: \
            [f'I={pth}' for pth in input['sorted_readgroup_bams']],
    threads: 4
    resources:
        mem_mb=8000,  # Merging sorted BAMs only takes <1GB of memory
        io_heavy=1
    log: 'logs/picard_merge_bams/{sample}.log'
    shell:
        "picard -Xmx{resources.mem_mb}m MergeSamFiles "
        "ASSUME_SORTED=false "
        "CREATE_INDEX=true "
        "MERGE_SEQUENCE_DICTIONARIES=false "
        "{params.input_bams_args} "
        "O={output.merged_bam} "
        "USE_THREADING=true "
        "SORT_ORDER=coordinate "
        "VALIDATION_STRINGENCY=STRICT "
        "TMP_DIR=$(mktemp -d --tmpdir tmp_picard.XXXXXXXXXXXXXXXX) "
        "2>{log} 1>&2"


rule picard_mark_dup_one_sample:
    """Picard MarkDuplicates on one sample."""
    output:
        bam='mark_dup_bam/{sample}.bam',
        bai='mark_dup_bam/{sample}.bai',
        metrics='mark_dup_bam/{sample}_marked_dup_metrics.txt'
    input:
        bam=rules.picard_merge_readgroup_bams_of_one_sample.output['merged_bam'],
        bai=rules.picard_merge_readgroup_bams_of_one_sample.output['merged_bai']
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 48000 + 8000 * (attempt - 1),
        io_heavy=1
    log: 'logs/picard_mark_dup/{sample}.log'
    shell:
        "picard -Xmx40000m MarkDuplicates "
        "CREATE_INDEX=true "
        "I={input.bam} O={output.bam} M={output.metrics} "
        "VALIDATION_STRINGENCY=STRICT "
        "TMP_DIR=$(mktemp -d --tmpdir tmp_picard.XXXXXXXXXXXXXXXX) "
        "2>{log} 1>&2"


rule copy_picard_bam_index:
    output: 'mark_dup_bam/{sample}.bam.bai'
    input: rules.picard_mark_dup_one_sample.output['bai']
    shell:
        "rsync -a {input} {output}"


def expand_to_all_samples(patterns):
    return {
        name: expand(pattern, sample=SAMPLES)
        for name, pattern in patterns.items()
    }


rule picard_mark_dup_all_samples:
    """Mark duplication for re-aligned BAMs of all samples."""
    input:
        **expand_to_all_samples({ \
            "sorted_bams": rules.picard_mark_dup_one_sample.output['bam'], \
            "sorted_bam_bais": rules.copy_picard_bam_index.output[0], \
        })


rule make_washu_output_manifest:
    """Generate the map of the custom aligned BAMs."""
    output:
        manifest='washu_dnaseq_alignment_summary.tsv'
    input: rules.picard_mark_dup_all_samples.input
    run:
        result_file_tpls = {
            'BAM': rules.picard_mark_dup_one_sample.output['bam'],
        }

        columns = [
            # This column is to be compaible to Matt's BAM map .dat file,
            # though its format is not the same.
            '# sample_name',
            'case', 'disease', 'experimental_strategy',
            'data_format', 'reference',
            'data_path', 'filesize',
            'UUID',  # also to be compatible to Matt's BAM map
        ]
        with open(output.manifest, 'w') as f:
            writer = csv.writer(f, dialect='excel-tab', lineterminator='\n')
            # Write column header
            writer.writerow(columns)

            # Write all generated output files
            for sample_name, sample_info in SAMPLE_INFO.items():
                for data_format, data_pth_fmt in result_file_tpls.items():
                    # Create sample_name column compatible to Matt's BAM map
                    new_sample_name = f"{sample_name}.hg38"

                    data_pth = Path(data_pth_fmt.format(sample=sample_name)).resolve()
                    file_size = data_pth.stat().st_size
                    writer.writerow([
                        new_sample_name,
                        sample_info.case_id,
                        sample_info.disease,
                        sample_info.experimental_strategy,
                        data_format, 'hg38',
                        str(data_pth), str(file_size),
                        "WUSTL-ADHOC-OUTPUT-NOT-TRACKED-BY-GDC"  # A fake UUID
                    ])
