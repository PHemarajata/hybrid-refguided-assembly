# Changelog

All notable changes to this project will be documented in this file.

---

## [v1.0.0] - Initial Release

### Overview
This is the **first public release** of the **Hybrid Reference-Guided Assembly Workflow**.  
The workflow was originally developed for **human monkeypox virus (hMPXV) Clade I**, which presents unique challenges due to its **long inverted terminal repeats (ITRs)**.  

### Features
- **Hybrid assembly support**: Illumina (short reads) + Oxford Nanopore (long reads).  
- **Reference-guided consensus**: using the Clade I reference genome `DQ011155.1`.  
- **ITR-aware masking**: includes a bundled 1-sided ITR BED (`references/itr-1sided.bed`).  
- **Automated pipeline**: implemented in Nextflow for reproducibility and resumability.  
- **QC reporting**: per-platform stats aggregated with MultiQC.  
- **Reproducible containerization**: via custom Docker image  
  `phemarajata614/viral-consensus:1.0.1`.  

### Inputs
- Illumina paired FASTQs (`--r1`, `--r2`)  
- ONT FASTQ or BAM (`--ont_fastq` or `--ont_bam`)  
- Reference genome FASTA (`--ref`)  
- ITR mask BED file (`--itr_bed`)  

### Outputs
- Aligned BAMs (Illumina, ONT, merged hybrid)  
- Variant calls (`calls.bcf`)  
- Mask file (`mask.bed`)  
- Consensus genome (`consensus.fa`)  
- QC metrics and a MultiQC report  

### Notes
- Best run with the **Docker profile** for reproducibility.  
- Initially designed for **hMPXV Clade I**, but may be adapted to other viral pathogens with similar genome structures.  

---
