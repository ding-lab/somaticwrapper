# SomaticWrapper (Tumor-Only) — compute1 v3.0

SomaticWrapper is a fully automated, modular pipeline for detecting **somatic variants from tumor-only WES/WGS** data on **GRCh38/HG38**.  
It runs on **LSF** (compute1) with per-chromosome parallelism and stage-to-stage job dependencies. Annotation uses **VEP 102**.

---

## Table of Contents
- [Features](#features)
- [Requirements](#requirements)
- [Install Third-Party Software](#install-third-party-software)
- [Input Layout](#input-layout)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Examples](#examples)
- [Tips & Troubleshooting](#tips--troubleshooting)
- [Acknowledgments](#acknowledgments)
- [Contact](#contact)
- [License](#license)

---

## Features
- Tumor-only variant calling on GRCh38/HG38
- **LSF-native** fan-out (per chromosome) and fan-in with `-w done(<JID>)`
- Mutect2 → merge/statistics/contamination → parsing/dbSNP filter
- **VEP 102** annotation & **MAF** generation per sample
- **Run-level report** after all samples finish

---

## Requirements
- **Platform:** LSF (compute1)
- **Docker images (used by the scripts):**
  - `scao/dailybox` (Mutect2, GATK, Picard, utilities)
  - `ensemblorg/ensembl-vep:release_102.0` (VEP 102 / vcf2maf)
- **References / Databases (HG38):**
  - GRCh38 FASTA (+ `.fai`)
  - gnomAD AF-only VCF
  - Panel of Normals (PoN)
  - Common biallelic sites list
  - VEP cache (release 102)
  - GTF (HG38)
- **Inputs:** one folder per sample under `--rdir`, each containing:
  - `<sample>.remDup.bam` (+ index)

> If your inputs are FASTQs, generate BAMs before running SomaticWrapper.

---

## Install Third-Party Software
Most tools are executed inside the Docker images above.  
If you maintain local installs, update the tool and reference paths at the top of `somaticwrapper.pl`.

---

## Input Layout
```
<run_dir>/
  SAMPLE_A/
    SAMPLE_A.remDup.bam
    SAMPLE_A.remDup.bam.bai
  SAMPLE_B/
    SAMPLE_B.remDup.bam
    SAMPLE_B.remDup.bam.bai
  ...
```

---

## Usage
From the directory where you cloned/downloaded **somaticwrapper**:

```bash
perl somaticwrapper.pl   --rdir <run_dir>   --ref  <path/to/GRCh38.fa>   --log  <log_dir>   --q    <queue>   --groupname <job_group_name>   --users <compute1_username>   --step <0..5>
```

### Arguments
| Flag | Description | Required | Default |
|---|---|---:|---|
| `--rdir` | Full path to run directory (contains per-sample subfolders) | ✅ | — |
| `--log` | Full path for logs; job scripts & LSF logs are written here | ✅ | — |
| `--ref` | GRCh38 FASTA (must have `.fai`) | ✅ | — |
| `--q` | LSF queue (e.g. `long`, `ding-lab`, `research-hpc`) |  | `long` |
| `--groupname` | Job group name (used as `/users/<users>/<groupname>`) | ✅ | — |
| `--users` | Compute1 username (used in group path) | ✅ | — |
| `--chr` | Reference has `chr` prefix (1=yes, 0=no) |  | `1` |
| `--sre` | Rerun mode (1=yes, 0=no) |  | `0` |
| `--exonic` | Exonic output flag for reporting (1=yes, 0=no) |  | `1` |
| `--step` | Stage to run (see below) | ✅ | — |

> **Job grouping:** the script submits with `-g /<users>/<groupname>`.  
> Example: `--users songcao --groupname SomaticWXS` → `-g /songcao/SomaticWXS`.

---

## Pipeline Steps
- **[0] Submit all steps with dependencies (recommended)**  
  For each sample:
  1. **j1**: Mutect2 per chromosome  
  2. **j2**: Gather VCFs + Merge stats + Read-orientation model + Contamination + FilterMutectCalls  
  3. **j3**: Parse & dbSNP filtering  
  4. **j4**: VEP 102 + vcf2maf → per-sample MAF  
  After **all j4** finish across samples → **j5** (run-level report).

- **[1]** Run Mutect2 (per-chromosome fan-out)  
- **[2]** Filter Mutect2 results (merge/statistics/contamination/filter)  
- **[3]** Parse Mutect2 results (coverage/VAF thresholds, dbSNP filter)  
- **[4]** Generate **per-sample** MAF (VEP 102 + vcf2maf)  
- **[5]** Generate **run-level** merged report (aggregates per-sample outputs)

> Step 6 (near-indel SNV removal / DNP annotate) has been **removed** in v3.0.

---

## Examples

**Run everything with dependencies (preferred):**
```bash
perl somaticwrapper.pl   --rdir /storage1/fs1/dinglab/Active/Projects/scao/gbm/GSAM   --log  /storage1/fs1/dinglab/Active/Projects/scao/gbm/GSAM.log   --ref  /storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa   --q long   --users songcao   --groupname SomaticWXS   --step 0
```

**Resume from filtering (step 2) only:**
```bash
perl somaticwrapper.pl ... --step 2
```

**Generate per-sample MAFs (step 4) only:**
```bash
perl somaticwrapper.pl ... --step 4
```

**Run the final merged report (step 5) only:**
```bash
perl somaticwrapper.pl ... --step 5
```

---

## Tips & Troubleshooting
- **Parallelism:** j1 shards by chromosome; downstream stages use `-w done(<JID>)` to enforce completion order.
- **Non-fast-forward pushes:** if you version the pipeline in Git and see a push rejection, rebase or merge remote changes before pushing.
- **Queues/Resources:** default queue is `long`; adapt the `bsub` resource strings if your environment differs.
- **Paths:** verify all reference/database paths at the top of `somaticwrapper.pl` match your environment.
- **VEP:** uses `ensemblorg/ensembl-vep:release_102.0` and a VEP 102 cache.

---

## Acknowledgments
- GATK / Mutect2, Picard, samtools
- Ensembl VEP (v102) & vcf2maf ecosystem
- gnomAD, COSMIC/DBSNP resources (where applicable)

---

## Contact
**Song Cao** — <scao@wustl.edu>

---

## License
Specify a license for this repository (e.g., MIT, BSD-3-Clause). Add a `LICENSE` file at the repo root.
