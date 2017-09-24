#!/bin/bash

# Create index files for reference:
# * faidx - required for ?
# * dict  - required for GATK
# * bwa index - required for alignment, skipping

# Usage: bash prepare_reference.sh FASTA
# Where FASTA is the full path to the reference file (FASTA format)

# Note that we assume FASTA ends in "fa"

PICARD="/usr/local/picard.jar"
SAMTOOLS="/usr/local/bin/samtools"
BWA="/usr/local/bwa"

function index_ref {
FA=$1

D=$(dirname $FA)
F=$(basename $FA)

CWD=`pwd`
cd $D

# DICT replaces "fa" suffix with "dict" - because GATK assumes F.fa has F.dict associated with it
DICT=${F%fa}dict

$SAMTOOLS faidx $F
#$BWA index $FA
java -jar $PICARD CreateSequenceDictionary R= $F O= $DICT

cd $CWD
}

index_ref $1
