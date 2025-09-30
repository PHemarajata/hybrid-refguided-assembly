#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Hybrid viral consensus (Illumina + ONT) with thread support
 * - ONT input as FASTQ(.gz) or BAM (BAM -> FASTQ)
 * - minimap2 mapping, bcftools calling/consensus
 * - FMT/DP disambiguation in filters
 * - Per-platform stats, MultiQC
 * - Fully resumable with per-step publishDir
 */

params.ref          = null
params.r1           = null
params.r2           = null
params.ont_fastq    = null
params.ont_bam      = null
params.sample       = 'sample'
params.min_dp       = 10
params.baseq        = 20
params.mapq         = 20
params.max_depth    = 100000
params.threads      = 120
params.outdir       = 'results'
params.itr_bed      = null
params.mask_lowcov  = true
params.multiqc_name   = "multiqc_report"
params.multiqc_outdir = "${params.outdir}/70_multiqc"


/********************
 * 00) Reference idx
 ********************/
process INDEX_REF {
  tag "${ref.baseName}"
  publishDir "${params.outdir}/00_ref", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path ref

  output:
  path "${ref}.fai"
  path ref

  script:
  """
  samtools faidx ${ref}
  """
}

/**************************
 * 10) Illumina mapping
 **************************/
process MAP_ILLUMINA {
  tag "${params.sample}"
  publishDir "${params.outdir}/10_map_illumina", mode: 'copy'
  cpus params.threads
  conda 'bioconda::minimap2=2.28 bioconda::samtools=1.20'

  input:
  path ref
  path fai
  path r1
  path r2

  output:
  path "illumina.sorted.bam"
  path "illumina.sorted.bam.bai"

  script:
  """
  minimap2 -t ${task.cpus} -ax sr ${ref} ${r1} ${r2} \
    | samtools sort -@ ${task.cpus} -o illumina.sorted.bam
  samtools index -@ ${task.cpus} illumina.sorted.bam
  """
}

/********************************
 * 20a) ONT BAM -> FASTQ (opt.)
 ********************************/
process ONT_BAM_TO_FASTQ {
  tag "${params.sample}"
  publishDir "${params.outdir}/20_map_ont", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path ontbam

  output:
  path "ont.extracted.fastq.gz"

  script:
  """
  samtools fastq -@ ${task.cpus} ${ontbam} | gzip -c > ont.extracted.fastq.gz
  """
}



/**************************
 * 20c) ONT mapping
 **************************/
process MAP_ONT {
  tag "${params.sample}"
  publishDir "${params.outdir}/20_map_ont", mode: 'copy'
  cpus params.threads
  conda 'bioconda::minimap2=2.28 bioconda::samtools=1.20'

  input:
  path ref
  path fai
  path ontfq

  output:
  path "ont.sorted.bam"
  path "ont.sorted.bam.bai"

  script:
  """
  minimap2 -t ${task.cpus} -ax map-ont ${ref} ${ontfq} \
    | samtools sort -@ ${task.cpus} -o ont.sorted.bam
  samtools index -@ ${task.cpus} ont.sorted.bam
  """
}

/************************
 * 30) Add read groups
 ************************/
process ADD_RG_ILLUMINA {
  tag "${params.sample}"
  publishDir "${params.outdir}/30_rg", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path bam

  output:
  path "illumina.rg.bam"
  path "illumina.rg.bam.bai"

  script:
  """
  samtools addreplacerg -@ ${task.cpus} -r ID:ILLUM -r PL:ILLUMINA -r SM:${params.sample} -o illumina.rg.bam ${bam}
  samtools index -@ ${task.cpus} illumina.rg.bam
  """
}

process ADD_RG_ONT {
  tag "${params.sample}"
  publishDir "${params.outdir}/30_rg", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path bam

  output:
  path "ont.rg.bam"
  path "ont.rg.bam.bai"

  script:
  """
  samtools addreplacerg -@ ${task.cpus} -r ID:ONT -r PL:ONT -r SM:${params.sample} -o ont.rg.bam ${bam}
  samtools index -@ ${task.cpus} ont.rg.bam
  """
}

/********************
 * 40) Merge BAM(s)
 ********************/
process MERGE_BAM {
  tag "${params.sample}"
  publishDir "${params.outdir}/40_merge", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path ref
  path fai
  path illum
  path ont, stageAs: 'ont.bam'

  output:
  path "hybrid.bam"
  path "hybrid.bam.bai"

  script:
  def ont_arg = ont.name != 'OPTIONAL_FILE' ? ont : ''
  """
  samtools merge -@ ${task.cpus} -o hybrid.bam ${illum} ${ont_arg}
  samtools index -@ ${task.cpus} hybrid.bam
  """
}

/***************************************
 * 50) Variants (unambiguous FMT/DP)
 ***************************************/
process CALL_VARIANTS {
  tag "${params.sample}"
  publishDir "${params.outdir}/50_variants", mode: 'copy'
  cpus params.threads
  conda 'bioconda::bcftools=1.20 bioconda::samtools=1.20'

  input:
  path ref
  path fai
  path bam

  output:
  path "calls.bcf"
  path "calls.bcf.csi"

  script:
  """
  bcftools mpileup -Ou \
    -f ${ref} \
    -a AD,ADF,ADR,DP \
    -Q ${params.baseq} -q ${params.mapq} \
    --max-depth ${params.max_depth} --threads ${task.cpus} \
    ${bam} \
  | bcftools call -mv --ploidy 1 -Ou \
  | bcftools filter -s LowCov -e 'FMT/DP<${params.min_dp}' -Ou \
  | bcftools norm -f ${ref} -Ob -o calls.bcf

  bcftools index --threads ${task.cpus} calls.bcf
  """
}

/****************************************
 * 55) Lowcov + optional ITR mask BEDs
 ****************************************/
process MAKE_MASKS {
  tag "${params.sample}"
  publishDir "${params.outdir}/50_variants", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20 bioconda::bedtools=2.31.1'

  input:
  path ref
  path bam
  path itr_bed, stageAs: 'itr.bed'

  output:
  path "lowcov.bed", optional: true
  path "mask.bed"

  script:
  def make_lowcov = params.mask_lowcov ? """
    samtools depth -@ ${task.cpus} -a -Q ${params.baseq} -q ${params.mapq} ${bam} \
      | awk '\$3<${params.min_dp} {print \$1"\\t"(\$2-1)"\\t"\$2}' > lowcov.bed
  """ : """
    : > lowcov.bed
  """
  def merge_cmd = itr_bed.name != 'OPTIONAL_FILE' ? """
    cat lowcov.bed ${itr_bed} | awk 'BEGIN{FS=OFS="\\t"} NF>=3 && \$2~/^[0-9]+\$/ && \$3~/^[0-9]+\$/ && \$3>\$2 {print \$1,\$2,\$3}' \
      | sort -k1,1 -k2,2n -k3,3n \
      | bedtools merge -i - > mask.bed
  """ : """
    awk 'BEGIN{FS=OFS="\\t"} NF==3 {print}' lowcov.bed \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - > mask.bed
  """
  """
  ${make_lowcov}
  ${merge_cmd}
  """
}

/*************************
 * 60) Consensus FASTA
 *************************/
process MAKE_CONSENSUS {
  tag "${params.sample}"
  publishDir "${params.outdir}/60_consensus", mode: 'copy'
  cpus params.threads
  conda 'bioconda::bcftools=1.20 bioconda::samtools=1.20'

  input:
  path ref
  path calls
  path mask

  output:
  path "consensus.fa"
  path "consensus.fa.fai"

  script:
  """
  # Ensure index is present/valid (belt & suspenders)
  bcftools index --threads ${task.cpus} -f ${calls}

  # Apply variants + masks to build haploid consensus
  bcftools consensus -f ${ref} -m ${mask} ${calls} > consensus_tmp.fa
  
  # Reheader the consensus FASTA with sample name as prefix
  # Get original reference name and prepend sample name
  original_name=\$(grep '^>' ${ref} | head -1 | sed 's/^>//')
  sed "s/^>.*/>${params.sample}_\$original_name/" consensus_tmp.fa > consensus.fa
  
  # Index the final consensus
  samtools faidx consensus.fa
  
  # Clean up temporary file
  rm consensus_tmp.fa
  """
}


/*************************
 * 52) bcftools stats
 *************************/
process BCFTOOLS_STATS {
  tag "${params.sample}"
  publishDir "${params.outdir}/50_variants", mode: 'copy'
  cpus params.threads
  conda 'bioconda::bcftools=1.20'

  input:
  path calls
  path ref

  output:
  path "bcftools.stats.txt"

  script:
  """
  bcftools stats --threads ${task.cpus} -F ${ref} ${calls} > bcftools.stats.txt
  """
}

/****************************************
 * 41) samtools stats (all reads)
 ****************************************/
process SAMTOOLS_STATS_ALL {
  tag "${params.sample}"
  publishDir "${params.outdir}/40_merge", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path ref
  path bam

  output:
  path "samtools.stats.all.txt"

  script:
  """
  samtools stats -@ ${task.cpus} -r ${ref} ${bam} > samtools.stats.all.txt
  """
}

/****************************************
 * 42) samtools stats/depth (ILLUM)
 ****************************************/
process SAMTOOLS_STATS_ILLUM {
  tag "${params.sample}"
  publishDir "${params.outdir}/40_merge", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path ref
  path bam

  output:
  path "samtools.stats.rg_ILLUM.txt"
  path "depth.rg_ILLUM.tsv"

  script:
  """
  samtools view -@ ${task.cpus} -u -r ILLUM ${bam} \
    | samtools stats -@ ${task.cpus} -r ${ref} - \
    > samtools.stats.rg_ILLUM.txt || true

  samtools view -@ ${task.cpus} -u -r ILLUM ${bam} \
    | samtools depth -@ ${task.cpus} -a -Q ${params.baseq} -q ${params.mapq} - \
    > depth.rg_ILLUM.tsv || true
  """
}

/****************************************
 * 43) samtools stats/depth (ONT)
 ****************************************/
process SAMTOOLS_STATS_ONT {
  tag "${params.sample}"
  publishDir "${params.outdir}/40_merge", mode: 'copy'
  cpus params.threads
  conda 'bioconda::samtools=1.20'

  input:
  path ref
  path bam

  output:
  path "samtools.stats.rg_ONT.txt"
  path "depth.rg_ONT.tsv"

  script:
  """
  samtools view -@ ${task.cpus} -u -r ONT ${bam} \
    | samtools stats -@ ${task.cpus} -r ${ref} - \
    > samtools.stats.rg_ONT.txt || true

  samtools view -@ ${task.cpus} -u -r ONT ${bam} \
    | samtools depth -@ ${task.cpus} -a -Q ${params.baseq} -q ${params.mapq} - \
    > depth.rg_ONT.tsv || true
  """
}

/****************
 * 70) MultiQC
 ****************/
// Default parameters for MultiQC
process MULTIQC {
  tag { "${params.multiqc_name}" }

  conda 'bioconda::multiqc=1.31'    // ← UNCOMMENT THIS

  publishDir params.multiqc_outdir, mode: 'copy', overwrite: true

  cpus 1
  memory '1 GB'
  time '30m'

  input:
  path samtools_stats_all
  path samtools_stats_rg_ILLUM
  path bcftools_stats
  path stats_ont_optional, stageAs: 'stats_ont.txt'

  output:
  path "*.html"
  path "*_data"

  script:
  def report_name = params.multiqc_name ?: 'multiqc_report'
  """
  set -euo pipefail

  mkdir -p multiqc_input

  # Copy files with standard names that MultiQC recognizes
  cp "${samtools_stats_all}" multiqc_input/samtools_stats.txt
  cp "${samtools_stats_rg_ILLUM}" multiqc_input/samtools_stats_illumina.txt  
  cp "${bcftools_stats}" multiqc_input/bcftools_stats.txt

  if [ "\$(basename "${stats_ont_optional}")" != "OPTIONAL_FILE" ]; then
    cp "${stats_ont_optional}" multiqc_input/samtools_stats_ont.txt
  fi

  echo "[MULTIQC] Input files:" >&2
  ls -la multiqc_input/ >&2
  
  echo "[MULTIQC] File sizes:" >&2
  wc -l multiqc_input/* >&2

  # Verify file content signatures
  echo "[MULTIQC] Checking samtools stats format:" >&2
  head -10 multiqc_input/samtools_stats.txt >&2
  
  echo "[MULTIQC] Checking bcftools stats format:" >&2  
  head -10 multiqc_input/bcftools_stats.txt >&2

  # Run MultiQC
  multiqc multiqc_input --force -n "${report_name}" --outdir . -v

  # Assert outputs exist
  test -f "${report_name}.html"
  test -d "${report_name}_data"
  """
}




workflow {
  // Parameter validation
  if (!params.ref) error "Please set --ref to your reference FASTA"
  if (!params.r1 || !params.r2) error "Please set --r1 and --r2 for Illumina data"
  if (!params.ont_fastq && !params.ont_bam) {
    log.warn "Proceeding with Illumina-only workflow (no ONT input)."
  }

  // Create input channels
  ch_ref_fa = Channel.fromPath(params.ref, checkIfExists: true)
  ch_r1 = Channel.fromPath(params.r1, checkIfExists: true)
  ch_r2 = Channel.fromPath(params.r2, checkIfExists: true)
  ch_ont_fq = params.ont_fastq ? Channel.fromPath(params.ont_fastq, checkIfExists: true) : Channel.empty()
  ch_ont_bam = params.ont_bam ? Channel.fromPath(params.ont_bam, checkIfExists: true) : Channel.empty()
  ch_itr_bed = params.itr_bed ? Channel.fromPath(params.itr_bed, checkIfExists: true) : Channel.value(file('OPTIONAL_FILE'))

  // Index reference
  INDEX_REF(ch_ref_fa)
  
  // Map Illumina reads
  MAP_ILLUMINA(INDEX_REF.out[1], INDEX_REF.out[0], ch_r1, ch_r2)
  
  // Handle ONT data if provided
  if (params.ont_bam) {
    ONT_BAM_TO_FASTQ(ch_ont_bam)
    ch_ont_fq_final = ONT_BAM_TO_FASTQ.out
  } else if (params.ont_fastq) {
    ch_ont_fq_final = ch_ont_fq
  } else {
    ch_ont_fq_final = Channel.empty()
  }
  
  // Map ONT reads if available
  if (params.ont_fastq || params.ont_bam) {
    MAP_ONT(INDEX_REF.out[1], INDEX_REF.out[0], ch_ont_fq_final)
    ch_ont_mapped = MAP_ONT.out[0]
  } else {
    ch_ont_mapped = Channel.empty()
  }
  
  // Add read groups
  ADD_RG_ILLUMINA(MAP_ILLUMINA.out[0])
  
  if (params.ont_fastq || params.ont_bam) {
    ADD_RG_ONT(ch_ont_mapped)
    ch_ont_rg = ADD_RG_ONT.out[0]
  } else {
    ch_ont_rg = Channel.value(file('OPTIONAL_FILE'))
  }
  
  // Merge BAMs
  MERGE_BAM(INDEX_REF.out[1], INDEX_REF.out[0], ADD_RG_ILLUMINA.out[0], ch_ont_rg)
  
  // Variant calling
  CALL_VARIANTS(INDEX_REF.out[1], INDEX_REF.out[0], MERGE_BAM.out[0])
  
  // Make masks
  MAKE_MASKS(INDEX_REF.out[1], MERGE_BAM.out[0], ch_itr_bed)
  
  // Make consensus
  MAKE_CONSENSUS(INDEX_REF.out[1], CALL_VARIANTS.out[0], MAKE_MASKS.out[1])
  
  // Statistics
  BCFTOOLS_STATS(CALL_VARIANTS.out[0], INDEX_REF.out[1])
  SAMTOOLS_STATS_ALL(INDEX_REF.out[1], MERGE_BAM.out[0])
  SAMTOOLS_STATS_ILLUM(INDEX_REF.out[1], MERGE_BAM.out[0])
  
  if (params.ont_fastq || params.ont_bam) {
    SAMTOOLS_STATS_ONT(INDEX_REF.out[1], MERGE_BAM.out[0])
    ch_ont_stats = SAMTOOLS_STATS_ONT.out[0]
  } else {
    ch_ont_stats = Channel.value(file('OPTIONAL_FILE'))
  }
  
  // MultiQC
  MULTIQC(
  SAMTOOLS_STATS_ALL.out,
  SAMTOOLS_STATS_ILLUM.out[0], 
  BCFTOOLS_STATS.out,              // ← Move this to 3rd position
  ch_ont_stats                     // ← Move this to 4th position 
  )
}
