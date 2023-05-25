import argparse
import os
import re
import logging
import subprocess
from pathlib import Path

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('sample', type=str,
    help='Sample id')

parser.add_argument('fq1', type=str,
    help='Fastq R1')

parser.add_argument('fq2', type=str,
    help='Fastq R2')

parser.add_argument('--flowcell', type=str, default='DUMMYFLOWCELL',
    help='flowcell used for sequencing')

parser.add_argument('--lane', type=str, default='DUMMYLANE',
    help='sequencing lane')

parser.add_argument('--index-sequencer', type=str, default='DUMMYSEQUENCER',
    help='index sequencer')

parser.add_argument('--library-preparation', type=str, default='DUMMYLIB',
    help='library prep id')

parser.add_argument('--known-sites', type=str,
    help='path to known sites file')

parser.add_argument('--reference', type=str,
    help='reference fp. Also needs to be in same directory as a .dict file')

parser.add_argument('--platform', type=str, default='ILLUMINA',
    help='platform. Default is ILLUMINA')

parser.add_argument('--out-prefix', type=str, default='output',
    help='output prefix for aligned and sorted bam file')

parser.add_argument('--cpu', type=int, default=1,
    help='cpus to use')


args = parser.parse_args()


def trimgalore(fq1, fq2, out_dir, cores=4, min_length=50):
    # $TRIMGALORE --phred33 --fastqc --cores 4 --length $MINLEN -q 20 -o $OUT --paired $FQ1 $FQ2 --path_to_cutadapt /diskmnt/Projects/Users/austins2/software/anaconda3/bin/cutadapt
    pieces = [
        'trim_galore',
        f'-phred33 --fastqc --cores {cores} -q 20',
        '--length', str(min_length),
        '-o', out_dir,
        '--paired', fq1, fq2
    ]
    return ' '.join(pieces)


def bwa_pe(sample, flowcell, lane, index_sequencer, library_preparation, platform,
            ref, fq1, fq2, out_sam, cpu=8):
    # $BWA mem -t 8 -M -R "@RG\tID:$NAME\tPL:illumina\tLB:$NAME\tPU:$NAME\tSM:$NAME" $REF_HUMAN $FQ1 $FQ2 | $SAMTOOLS view -Shb -o $OUT/$NAME.human.bam -
    id = f'{flowcell}.{lane}'
    pl = platform
    lb = f'{sample}.{library_preparation}'
    pu = f'{flowcell}.{lane}.{index_sequencer}'
    sm = sample

    pieces = [
        'bwa mem',
        '-t', str(cpu),
        '-M',
        '-R', f'"@RG\\tID:{id}\\tPL:{pl}\\tLB:{lb}\\tPU:{pu}\\tSM:{sm}"',
        '-o', out_sam,
        ref, fq1, fq2
    ]
    return ' '.join(pieces)


def sam_to_bam(sam_fp, bam_fp):
    pieces = [
        'samtools view -hb',
        '-o', bam_fp,
        sam_fp
    ]
    return ' '.join(pieces)


def sort_and_index(input_bam, output_bam):
    # $JAVA -Xmx16G -jar $PICARD SortSam \
    #    CREATE_INDEX=true \
    #    I=$OUT/$NAME.human.bam \
    #    O=$OUT/$NAME.human.sorted.bam \
    #    SORT_ORDER=coordinate \
    #    VALIDATION_STRINGENCY=STRICT
    pieces = [
        'picard SortSam CREATE_INDEX=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT',
        f'I={input_bam}',
        f'O={output_bam}',
    ]
    return ' '.join(pieces)


def remove_duplicates(input_bam, output_bam, metrics_fp):
    # # remove-duplication
    # echo "[INFO] s3: picard - markduplicates" >&2
    # $JAVA -Xmx16G -jar $PICARD MarkDuplicates \
    #    I=$OUT/$NAME.human.sorted.bam \
    #    O=$OUT/$NAME.human.remDup.bam \
    #    REMOVE_DUPLICATES=true \
    #    M=$OUT/$NAME.human.remDup.metrics.txt
    pieces = [
        'picard MarkDuplicates REMOVE_DUPLICATES=false',
        f'I={input_bam}',
        f'O={output_bam}',
        f'M={metrics_fp}',
    ]
    return ' '.join(pieces)


def base_recalibrator(input_bam, known_sites, ref, bsqr_file):
    pieces = [
        'gatk BaseRecalibrator',
        '--input', input_bam,
        '--known-sites', known_sites,
        '--reference', ref,
        '--output', bsqr_file
    ]
    return ' '.join(pieces)


def apply_base_recalibrator(input_bam, bsqr_file, output_fp):
    pieces = [
        'gatk ApplyBQSR',
        '--input', input_bam,
        '--bqsr-recal-file', bsqr_file,
        '--emit-original-quals', 'true',
        '--output', output_fp
    ]
    return ' '.join(pieces)


def index_bam(input_bam):
    return f'samtools index {input_bam}'


def run_align_dnaseq(fq1, fq2, reference, known_sites, sample, flowcell, lane, index_sequencer, library_preparation, platform, output_prefix, cpu):
    logging.info('running trim galore')
    trim_galore_out_dir = 'trim_galore_outputs'
    cmd = trimgalore(fq1, fq2, trim_galore_out_dir)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    fq1_root = re.sub(r'^(.*).((fq)|(fastq))(.gz)?', r'\1', fq1.split('/')[-1])
    fq2_root = re.sub(r'^(.*).((fq)|(fastq))(.gz)?', r'\1', fq2.split('/')[-1])
    trimmed_fq1 = os.path.join(trim_galore_out_dir, f'{fq1_root}_val_1.fq.gz')
    trimmed_fq2 = os.path.join(trim_galore_out_dir, f'{fq2_root}_val_2.fq.gz')

    logging.info('aligning with bwa')
    intermediate_dir = 'intermediates'
    Path(intermediate_dir).mkdir(exist_ok=True, parents=True)
    out_sam = os.path.join(intermediate_dir, 'bwa_out.sam')
    cmd = bwa_pe(sample, flowcell, lane, index_sequencer, library_preparation,
                 platform, reference, trimmed_fq1, trimmed_fq2, out_sam, cpu=8)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    logging.info('converting sam to bam')
    bwa_bam = os.path.join(intermediate_dir, 'bwa_out.bam')
    cmd = sam_to_bam(out_sam, bwa_bam)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    logging.info('sorting and indexing bam')
    sorted_bam = os.path.join(intermediate_dir, 'bwa_out.sorted.bam')
    cmd = sort_and_index(bwa_bam, sorted_bam)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    logging.info('removing duplicates')
    dedup_bam = os.path.join(intermediate_dir, 'bwa_out.sorted.dedup.bam')
    dedup_metrics = os.path.join(intermediate_dir, 'dedup_metrics.txt')
    cmd = remove_duplicates(sorted_bam, dedup_bam, dedup_metrics)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    logging.info('modeling bsqr')
    bsqr_file = os.path.join(intermediate_dir, 'bsqr_recal_file.table')
    cmd = base_recalibrator(dedup_bam, known_sites, reference, bsqr_file)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    logging.info('applying bsqr')
    bsqr_bam = os.path.join(intermediate_dir, 'final.bam')
    cmd = apply_base_recalibrator(dedup_bam, bsqr_file, bsqr_bam)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    logging.info('sorting and indexing bam')
    output_bam = f'{output_prefix}.bam'
    cmd = sort_and_index(bsqr_bam, output_bam)
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info(output)

    # picard outputs index as just output.bai not output.bam.bai so renaming
    os.rename(f'{output_prefix}.bai', f'{output_prefix}.bam.bai')

    logging.info('cleaning up large intermediates')
    for fp in [out_sam, bwa_bam, sorted_bam, dedup_bam, bsqr_bam]:
        logging.info(f'removing {fp}')
        os.remove(fp)


def main():
    run_align_dnaseq(
        args.fq1, args.fq2, args.reference, args.known_sites, args.sample,
        args.flowcell, args.lane, args.index_sequencer,
        args.library_preparation, args.platform, args.out_prefix, args.cpu)


if __name__ == '__main__':
    main()
