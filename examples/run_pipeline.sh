#!/usr/bin/env bash
set -euo pipefail

OUTDIR="results"
SAMPLES="config/samples.csv"
THREADS=16

# 1. Raw read QC
trebouxia-pipeline --outdir "$OUTDIR" qc --samples "$SAMPLES" --threads "$THREADS"

# 2. Adapter and quality trimming
trebouxia-pipeline --outdir "$OUTDIR" trim --samples "$SAMPLES" --threads "$THREADS"

# 3. SPAdes assemblies
trebouxia-pipeline --outdir "$OUTDIR" assemble --samples "$SAMPLES" --threads 32

# 4. Assembly quality checks
trebouxia-pipeline --outdir "$OUTDIR" assembly-qc --samples "$SAMPLES" --threads "$THREADS" --busco-lineage chlorophyta_odb10

# 5. Whole-genome comparison
trebouxia-pipeline --outdir "$OUTDIR" compare --samples "$SAMPLES" --threads "$THREADS"

# 6. Read-backed variants: ASV2 reads mapped to ASV1 assembly
trebouxia-pipeline --outdir "$OUTDIR" variants \
  --reference "$OUTDIR/02_assembly/ASV1/scaffolds.fasta" \
  --query-sample ASV2 \
  --samples "$SAMPLES" \
  --threads "$THREADS" \
  --min-qual 30 \
  --min-depth 10

# 7. SnpEff annotation, assuming a local database called Tlynniae_ASV1
trebouxia-pipeline --outdir "$OUTDIR" snpeff \
  --vcf "$OUTDIR/06_variants/ASV2_vs_ASV1.filtered.vcf.gz" \
  --snpeff-db Tlynniae_ASV1

# 8. Summarise HIGH and MODERATE impact variants
trebouxia-pipeline --outdir "$OUTDIR" summarise \
  --ann-vcf "$OUTDIR/06_variants/snpeff/variants.ann.vcf" \
  --impact HIGH,MODERATE
