#Trebouxia_Genome_Comparition
# Trebouxia strain comparative genomics pipeline

Python CLI pipeline for comparing two closely related *Trebouxia lynniae* strains at genome and variant-impact level.

The workflow follows the genomics logic used in the manuscript: Illumina read QC/trimming, SPAdes assemblies, QUAST/BUSCO assembly checks, FastANI and MUMmer whole-genome comparison, SnpEff variant-effect annotation, and integration with EggNOG-mapper and InterProScan functional annotations.

## Repository layout

```text
trebouxia_genomics_pipeline/
├── trebouxia_pipeline/
│   ├── cli.py
│   ├── merge_annotations.py
│   └── plots.py
├── config/
│   ├── samples.csv
│   └── environment.yml
├── examples/
│   └── run_pipeline.sh
├── tests/
│   └── test_snpeff_parser.py
├── pyproject.toml
├── .gitignore
├── LICENSE
└── README.md
```

## Installation

```bash
git clone https://github.com/YOUR_USER/trebouxia_genomics_pipeline.git
cd trebouxia_genomics_pipeline
mamba env create -f config/environment.yml
mamba activate trebouxia-genomics
pip install -e .
```

Some heavy tools are best installed separately depending on your cluster: MAKER, InterProScan, EggNOG-mapper databases and a local SnpEff database built from the reference annotation.

## Input

Edit `config/samples.csv`:

```csv
sample,r1,r2
ASV1,/path/ASV1_R1.fastq.gz,/path/ASV1_R2.fastq.gz
ASV2,/path/ASV2_R1.fastq.gz,/path/ASV2_R2.fastq.gz
```

## Minimal run

```bash
trebouxia-pipeline --outdir results qc --samples config/samples.csv --threads 16
trebouxia-pipeline --outdir results assemble --samples config/samples.csv --threads 32
trebouxia-pipeline --outdir results compare --samples config/samples.csv --threads 32
```

## Variant-effect workflow

SnpEff requires a VCF and a custom database matching your reference genome annotation. A conservative option is to map ASV2 reads to ASV1, call variants with `bcftools`, and annotate them against the ASV1 gene models.

```bash
trebouxia-pipeline --outdir results snpeff \
  --vcf results/06_variants/ASV2_vs_ASV1.filtered.vcf.gz \
  --snpeff-db Tlynniae_ASV1

trebouxia-pipeline --outdir results summarise \
  --ann-vcf results/06_variants/snpeff/variants.ann.vcf \
  --impact HIGH,MODERATE
```

## Functional annotation merge

```bash
python -m trebouxia_pipeline.merge_annotations \
  --gene-counts results/07_functional_summary/variant_counts_by_gene.tsv \
  --eggnog annotations/ASV1.emapper.annotations \
  --interpro annotations/ASV1.interproscan.tsv \
  --output results/07_functional_summary/variant_genes_functional_annotation.tsv
```

## Quick plots

```bash
python -m trebouxia_pipeline.plots \
  --input results/07_functional_summary/variant_effects_long.tsv \
  --column effect \
  --output results/07_functional_summary/snpEff_effects_barplot.png
```

## Recommended strict genomics additions before publication

1. Perform read-backed variant calling rather than relying only on assembly-to-assembly SNPs.
2. Mask repeats and low-complexity regions before interpreting private/exclusive genes.
3. Report callable genome size and coverage distribution for both strains.
4. Separate nuclear, chloroplast and mitochondrial comparisons.
5. For “exclusive genes”, require read coverage support and reciprocal protein searches.
6. Report whether non-synonymous variants are homozygous/fixed and avoid overinterpreting “high impact” calls without manual inspection.
7. Link key genes to phenotypes only as candidates unless validated by RNA-seq, proteomics or allele-specific assays.

## Citation

Please cite the manuscript and the underlying tools used by each step.
