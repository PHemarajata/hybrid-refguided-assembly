# Hybrid Viral Consensus Pipeline (Illumina + ONT)  

A reproducible, modular Nextflow (DSL2) workflow to generate a **consensus viral genome** from **Illumina short reads** and optional **ONT (long) reads or BAM**, including coverage masking, terminal repeat handling, and QC reporting.

---

## 🚀 Purpose

- Generate a high‑confidence consensus sequence for a viral sample (e.g., mpox) by combining Illumina + ONT data.  
- Mask low‑coverage and repeat (ITR) regions to avoid misleading confident calls in ambiguous zones.  
- Provide QC metrics and summaries via **MultiQC**, per-platform depth stats, variant counts, etc.  
- Be fully resumable, container‑aware (conda, Docker, Singularity), and thread‑efficient.

---

## 🔍 Features

- Accepts ONT input as **FASTQ(.gz)** or **BAM** (converts BAM → FASTQ internally).  
- Maps Illumina (`-x sr`) and ONT (`-x map-ont`) reads with **minimap2**, sorts and indexes BAMs.  
- Adds Read Groups (RG: ILLUM, ONT) for platform-aware stats and filtering.  
- Merges alignments into a **hybrid BAM**.  
- Calls variants with **bcftools** (haploid model) and **disambiguated DP filtering** (`FMT/DP < min_dp`).  
- Builds a **mask BED** combining **low-coverage** regions + optional ITR BED merged  
- Generates **consensus FASTA**, replacing masked regions with `N`.  
- Computes platform-specific depth & stats (samtools), variant stats (bcftools), and aggregates via **MultiQC**.  
- Supports multiple execution profiles: **conda**, **Docker**, **Singularity**.  
- Fully resume-able: each logical step is isolated with input/output channels and caching.

---

## 📦 Repository Structure

```
.
├── main.nf                ← Nextflow pipeline
├── nextflow.config        ← Configuration (threads, containers, profiles)
├── README.md              ← This document
└── workflows/             ← (optional) example input / test scenarios
```

---

## 🧩 Workflow Overview

1. **Index Reference** (faidx)  
2. **Map Illumina Reads** → sorted BAM  
3. **Map ONT Reads** (if present; BAM → FASTQ → align) → sorted BAM  
4. **Add Read Groups** to Illumina & ONT BAMs  
5. **Merge** BAMs into a hybrid BAM  
6. **Call Variants** using bcftools, using `FMT/DP`-based filtering  
7. **Make Masks**: low-coverage BED + optional ITR BED merged  
8. **Generate Consensus**: patch variants + mask = consensus FASTA  
9. **Compute Stats**: samtools / bcftools per‑RG and combined  
10. **MultiQC**: aggregate QC & stats into HTML report  

---

## 📥 Usage

```bash
nextflow run main.nf   --ref /path/to/ref.fa   --r1 /path/to/sample_R1.fastq.gz   --r2 /path/to/sample_R2.fastq.gz   [--ont_fastq /path/to/ont.fastq.gz | --ont_bam /path/to/ont.bam]   [--itr_bed /path/to/itr_mask.bed]   --sample sample_name   --threads 16   --min_dp 10   --outdir results
```

- Use `-profile docker` or `-profile singularity` to run within containers.  
- Use `-resume` for incremental workflow runs.  
- Example (Illumina-only):
  
  ```bash
  nextflow run main.nf     --ref ref.fa     --r1 R1.gz     --r2 R2.gz     --sample mysample     --outdir results_illumina
  ```

---

## 🗄️ Inputs

| Name         | Description |
|--------------|-------------|
| `--ref`       | Reference FASTA (e.g. DQ011155.1) |
| `--r1` / `--r2` | Illumina paired FASTQ (gzip OK) |
| `--ont_fastq` | ONT long reads FASTQ (.gz) |
| `--ont_bam`   | ONT-aligned BAM (alternative to FASTQ) |
| `--itr_bed`   | BED file of terminal repeat region(s) to be masked |

---

## 📦 Outputs (under `${outdir}`)

```
outdir/
  00_ref/           ← reference + .fai  
  10_map_illumina/  ← illumina.sorted.bam(.bai)  
  20_map_ont/       ← ont.sorted.bam(.bai) or extracted FASTQ  
  30_rg/            ← illumina.rg.bam(.bai), ont.rg.bam(.bai)  
  40_merge/         ← hybrid.bam(.bai), per-RG stats, depth files  
  50_variants/      ← calls.bcf(.csi), mask.bed, lowcov.bed  
  60_consensus/     ← consensus.fa(.fai)  
  70_multiqc/       ← multiqc_report.html + multiqc_data  
```

---

## 📝 Notes & Caveats

- **Masking is essential**. Without a `mask.bed`, consensus will fill zero-coverage regions with reference, misleading downstream analysis.  
- The VCF only captures **small variants** — large structural changes may not be included.  
- ITRs (terminal repeats) can complicate mapping; mask or handle them explicitly.  
- Tools like `bcftools consensus` require **indexed VCF/BCF**; we stage and/or (re)index to ensure that.  
- MultiQC only writes a report if it detects at least one supported file in its input directory. Ensure your stats files are staged correctly.  
- If ONT data is absent, the pipeline gracefully skips ONT steps and runs Illumina-only mode.  
- The workflow assumes **single-sample, haploid viral genome**. Mixed infections or contamination may yield odd consensus calls.  
- Be mindful of resource usage: `--threads` applies to mapping, indexing, and variant calling. Adjust per your machine.

---

## 🧪 Testing & Validation

- You can provide a small test dataset (e.g. subset of reads + known reference) in `workflows/test/` and run:

  ```bash
  nextflow run main.nf -profile test --ref test_ref.fa --r1 test_R1.gz --r2 test_R2.gz
  ```

- Use `-with-dag flowchart.svg` to visualize step dependencies.  
- Inspect `work/<hash>/` for `.command.sh` and logs when a step fails.

---

## ✅ Summary

This pipeline gives you a robust, reproducible way to go from raw Illumina + ONT (or ONT BAM) to a high-confidence, masked viral consensus sequence, plus QC summaries and platform-aware depth metrics. It is portable, resumable, and built for typical viral genomics conditions like mpox (ITRs, low coverage, mixed data types).
