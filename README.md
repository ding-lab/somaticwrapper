# SomaticWrapper v3.0.1 (Compute1)

### Automated Somatic Variant Calling Pipeline (HG38)

**SomaticWrapper** is a fully automated and modular pipeline for detecting somatic variants from paired tumor–normal WGS/WXS data on the LSF compute1 cluster (WashU).  
It integrates multiple industry-standard variant callers — **Strelka2**, **VarScan2**, **Mutect1**, and **Pindel** — and produces comprehensive, annotated mutation calls in MAF format.

---

## 🔬 Overview

- **SNV calls**: intersection of *2 out of 3* callers — Strelka2, Mutect1, VarScan2  
- **Indel calls**: intersection of *2 out of 3* callers — Strelka2, VarScan2, Pindel  
- **Reference genome**: Human GRCh38 (HG38)  
- **Scheduler**: LSF (supports job dependencies and groups)

Final output files:
- `dnp.annotated.maf` → all variants  
- `dnp.annotated.coding.maf` → coding variants only  

---

## 🚀 Improvements (v3.0.1)

1. Added sample matching checking (step 15)

---

## ⚙️ Environment Setup (Compute1)

Before running, update your `~/.bashrc` to include the necessary environment variables:

```bash
export PATH=/storage1/fs1/songcao/Active/Software/anaconda3/bin:$PATH
export STORAGE1=/storage1/fs1/songcao/Active
export STORAGE2=/storage1/fs1/dinglab/Active
export STORAGE3=/storage1/fs1/m.wyczalkowski/Active
export LSF_DOCKER_VOLUMES="$STORAGE1:$STORAGE1 $STORAGE2:$STORAGE2 $STORAGE3:$STORAGE3"
```

Then activate:
```bash
source ~/.bashrc
```

---

## 🧩 Usage

### Step 1. Download or clone this repository
```bash
git clone https://github.com/YourGitRepo/somaticwrapper.git
cd somaticwrapper
```

### Step 2. Prepare your run and log directories
Example:
```bash
mkdir -p /storage1/fs1/songcao/Active/Projects/somatic/example_run_somatic_2025
mkdir -p /storage1/fs1/songcao/Active/Projects/somatic/example_run_somatic_2025/log
```

### Step 3. Run the pipeline

#### Option A – Run all steps automatically (recommended)
Use **`--step 0`** to run Steps 1–14 sequentially with built-in job dependencies:
```bash
perl somaticwrapper.pl   --step 0   --rdir /storage1/fs1/songcao/Active/Projects/somatic/example_run_somatic_2025   --log  /storage1/fs1/songcao/Active/Projects/somatic/example_run_somatic_2025/log   --ref /storage1/fs1/dinglab/Active/Projects/scao/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa   --smg /storage1/fs1/dinglab/Active/Projects/scao/Database/SMG/smg_list.txt   --groupname example_run_somatic_2025   --users scao   --wgs 0   --srg 1   --sre 0   --exonic 1   --q long   --mincovt 14 --mincovn 8 --minvaf 0.05 --maxindsize 100
```

#### Option B – Run a specific step manually
```bash
perl somaticwrapper.pl --step 5 --rdir <run_dir> --log <log_dir> ...
```

---

## 🔢 Step Reference

| Step | Description |
|------|--------------|
| 0 | **Submit steps (1–11) automatically with dependencies** |
| 1 | Run Strelka2 |
| 2 | Run VarScan2 |
| 3 | Run Pindel |
| 4 | Run Mutect1 |
| 5 | Parse Mutect results |
| 6 | Parse Strelka2 results |
| 7 | Parse VarScan2 results |
| 8 | Parse Pindel results |
| 9 | QC VCF files |
| 10 | Merge VCF files |
| 11 | Generate MAF files |
| 12 | Merge run-level MAF |
| 13 | DNP annotation |
| 14 | Clean unnecessary intermediate files |
| 22 | **Submit steps (2–11) automatically with dependencies** |
| 23 | **Submit steps (3–11) automatically with dependencies** |
| 24 | **Submit steps (4–11) automatically with dependencies** |
| 25 | **Submit steps (5–11) automatically with dependencies** |
| 26 | **Submit steps (6–11) automatically with dependencies** |
| 27 | **Submit steps (7–11) automatically with dependencies** |
| 28 | **Submit steps (8–11) automatically with dependencies** |
| 29 | **Submit steps (9–11) automatically with dependencies** |
| 30 | **Submit steps (10–11) automatically with dependencies** |
---

## ⚙️ Key Parameters

| Parameter | Description |
|------------|-------------|
| `--rdir` | Full path to run directory containing per-sample folders |
| `--log` | Path for log output (usually parent of rdir) |
| `--srg` | BAM has read groups (1 = yes, 0 = no) |
| `--sre` | Rerun and overwrite results (1 = yes, 0 = no) |
| `--wgs` | 1 = WGS, 0 = WXS |
| `--groupname` | Job group name |
| `--users` | LSF user account (used in job group path) |
| `--ref` | HG38 reference FASTA |
| `--smg` | SMG gene list file |
| `--q` | LSF queue (`research-hpc`, `ding-lab`, or `long`) |
| `--mincovt` | Minimum tumor coverage (≥ 14) |
| `--mincovn` | Minimum normal coverage (≥ 8) |
| `--minvaf` | Minimum variant allele frequency (≥ 0.05) |
| `--maxindsize` | Maximum indel size (≤ 100) |
| `--exonic` | Output exonic region (1 = yes, 0 = no) |

---

## 🧾 Example Output Files

```
run_dir/
├── <sample_name>/
│   ├── strelka/
│   ├── varscan/
│   ├── pindel/
│   ├── mutect1/
│   ├── merged.withmutect.vcf
│   ├── <sample>.withmutect.maf
│   └── <sample>.dnp.annotated.maf
└── log/
    ├── LSF_DIR_SOMATIC/
    └── tmpsomatic/
```

---

## 👤 Contact

**Author:** Song Cao  
**Email:** [scao@wustl.edu](mailto:scao@wustl.edu)  
Washington University in St. Louis 

---
