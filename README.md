# Hybrid Reference-Guided Assembly Workflow

## 🎯 Purpose
This workflow performs **reference-guided consensus genome assembly** for viral datasets (e.g., **mpox Clade I**), supporting both **Illumina short reads** and **Oxford Nanopore long reads**.  

Key goals:
- Generate a **haploid consensus genome**.
- **Mask low coverage and repetitive ITR regions**.
- Provide **QC summaries** for each platform.
- Ensure **reproducibility** with containerization.
- Remain **resumable** and **portable** across environments.

---

## ✨ Features
- Supports **Illumina (FASTQ R1/R2)** + **ONT (FASTQ or BAM)** simultaneously.
- Flexible read mapping with **minimap2** (`-x sr` and `-x map-ont` modes).
- Read group tagging for platform separation (`ILLUM`, `ONT`).
- Variant calling with **bcftools mpileup/call**, using disambiguated `DP` handling.
- Mask creation:
  - **Low coverage** sites (configurable cutoff, default DP<10).
  - **ITR regions** (reference BED provided).
- Consensus genome generation with **bcftools consensus**.
- Per-platform stats with **samtools** and aggregation in **MultiQC**.
- Fully resumable via Nextflow `-resume`.
- Reproducible execution with **Docker container**:  
  `phemarajata614/viral-consensus:1.0.1`.

---

## 📂 Repository Structure
```
├── main.nf             # Main Nextflow pipeline
├── nextflow.config     # Configuration (Docker, threads, defaults)
├── README.md           # This document
└── references/         # Bundled reference & ITR mask
    ├── DQ011155.1.fasta
    └── itr-1sided.bed
```

---

## 🔄 Workflow Overview
1. **Reference indexing** (`samtools faidx`)
2. **Read mapping**
   - Illumina → minimap2 (`-x sr`)
   - ONT → minimap2 (`-x map-ont`)
3. **Read group tagging** (`samtools addreplacerg`)
4. **Merge BAMs** into hybrid alignment
5. **Variant calling** (`bcftools mpileup/call`)
6. **Mask generation** (low coverage + ITR BED)
7. **Consensus building** (`bcftools consensus`)
8. **Stats collection** (`samtools stats`, `depth`)
9. **MultiQC report** summarizing results

---

## 🚀 Quick Start

**Prereqs:** Nextflow ≥ 23.x and Docker installed.

**Container:**  
`phemarajata614/viral-consensus:1.0.1`

**Bundled references:**  
- `references/DQ011155.1.fasta` → Mpox Clade I reference  
- `references/itr-1sided.bed` → 1-sided ITR mask  

### Example (Hybrid Illumina + ONT)
```bash
nextflow run main.nf -profile laptop \
  --ref references/DQ011155.1.fasta \
  --r1  /path/to/R1.fastq.gz \
  --r2  /path/to/R2.fastq.gz \
  --ont_fastq /path/to/ont.fastq.gz \
  --itr_bed references/itr-1sided.bed \
  --sample my_sample \
  --outdir data/my_sample/results
```

### Example (Illumina-only)
```bash
nextflow run main.nf -profile laptop \
  --ref references/DQ011155.1.fasta \
  --r1  /path/to/R1.fastq.gz \
  --r2  /path/to/R2.fastq.gz \
  --itr_bed references/itr-1sided.bed \
  --sample my_sample \
  --outdir data/my_sample/results
```

---

## 🔧 Hardware Profiles

Select the profile that matches your machine with `-profile <name>`. Thread count and memory limits are set automatically — no need to pass `--threads`.

| Profile | Target machine | CPUs | Memory | maxForks |
| ------- | ------------- | ---- | ------ | -------- |
| `laptop` | 22-core workstation / laptop (62 GB RAM) | 18 | 52 GB | 2 |
| `a100` | A100 workstation (128-core, 512 GB RAM) | 120 | 480 GB | 1 |
| `docker` | Generic / CI (override with `--threads`) | 16 | 16 GB | 4 |

```bash
# Laptop / local workstation
nextflow run main.nf -profile laptop ...

# A100 high-memory workstation
nextflow run main.nf -profile a100 ...
```

---

## 🧾 Inputs

| Param         | Required | Notes |
|---------------|----------|-------|
| `--ref`       | ✅ | Reference FASTA; Mpox Clade I provided. |
| `--r1`, `--r2`| ✅ | Illumina paired FASTQs (`.gz` OK). |
| `--ont_fastq` | optional | ONT long-read FASTQ. |
| `--ont_bam`   | optional | ONT BAM; auto-converted to FASTQ. |
| `--itr_bed`   | optional | Mask BED; ITR BED provided. |
| `--sample`    | ✅ | Sample name (used in RG tags + outputs). |
| `--threads`   | ✅ | CPU threads (propagated to key tools). |
| `--outdir`    | ✅ | Output root directory. |

---

## 📦 Outputs

- `00_ref/` → reference + index  
- `10_map_illumina/` → Illumina BAMs  
- `20_map_ont/` → ONT BAM/FASTQ  
- `30_rg/` → read-group annotated BAMs  
- `40_merge/` → merged hybrid BAM + stats  
- `50_variants/` → calls.bcf, mask.bed  
- `60_consensus/` → consensus.fa + index  
- `70_multiqc/` → MultiQC report  

---

## ⚠️ Notes & Caveats
- **Masking is essential**: without BED, consensus may double-count ITRs.  
- **VCF-based consensus**: captures SNPs/indels but not structural events.  
- **MultiQC** will only produce reports if expected stats files exist.  
- **Single-sample assumption**: designed for haploid viral genomes.  

---

## 🧪 Testing & Validation
- Try with bundled **mpox Clade I reference** + **test datasets**.  
- Generate workflow DAG:  
  ```bash
  nextflow run main.nf -profile docker -with-dag flowchart.png
  ```  
- Inspect task scripts:  
  ```bash
  less work/<hash>/.command.sh
  ```  

---

## ✅ Summary
This workflow provides a **robust, reproducible, containerized method** for generating viral consensus genomes from hybrid Illumina + ONT data. It is especially suited to **mpox Clade I** with **ITR complexity**, but the design is generalizable.  
