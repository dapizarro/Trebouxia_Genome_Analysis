# Trebouxia_Genome_Comparison

> Trebouxia strain comparative genomics pipeline

Python CLI pipeline for comparing two closely related *Trebouxia lynniae* strains at genome and variant-impact level.

The workflow follows a strict comparative genomics logic: Illumina read QC/trimming, SPAdes assemblies, QUAST/BUSCO assembly checks, FastANI and MUMmer whole-genome comparison, read-backed variant calling, SnpEff variant-effect annotation, and integration with EggNOG-mapper and InterProScan functional annotations.

---

## Repository layout

```text
trebouxia_genomics_pipeline/
├── trebouxia_pipeline/
│   ├── __init__.py
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

---

## Installation

```bash
git clone https://github.com/dapizarro/trebouxia_genomics_pipeline.git
cd trebouxia_genomics_pipeline
mamba env create -f config/environment.yml
mamba activate trebouxia-genomics
pip install -e .
```

Some heavy tools are best installed separately depending on your cluster: MAKER, InterProScan, EggNOG-mapper databases and a local SnpEff database built from the reference annotation.

---

## Input

Edit `config/samples.csv`:

```csv
sample,r1,r2
ASV1,/path/ASV1_R1.fastq.gz,/path/ASV1_R2.fastq.gz
ASV2,/path/ASV2_R1.fastq.gz,/path/ASV2_R2.fastq.gz
```

For genome-to-genome comparisons, provide assembled genomes after the `assemble` step or manually in:

```text
results/02_assembly/ASV1/scaffolds.fasta
results/02_assembly/ASV2/scaffolds.fasta
```

---

## Minimal run

```bash
trebouxia-pipeline --outdir results qc --samples config/samples.csv --threads 16
trebouxia-pipeline --outdir results assemble --samples config/samples.csv --threads 32
trebouxia-pipeline --outdir results compare --samples config/samples.csv --threads 32
```

---

## Read-backed variant workflow

A conservative option is to map ASV2 reads to ASV1, call variants with `bcftools`, filter them, and annotate them against the ASV1 gene models using a local SnpEff database.

```bash
trebouxia-pipeline --outdir results variants \
  --reference results/02_assembly/ASV1/scaffolds.fasta \
  --query-sample ASV2 \
  --samples config/samples.csv \
  --threads 32
```

---

## SnpEff variant-effect workflow

SnpEff requires a VCF and a custom database matching your reference genome annotation.

```bash
trebouxia-pipeline --outdir results snpeff \
  --vcf results/06_variants/ASV2_vs_ASV1.filtered.vcf.gz \
  --snpeff-db Tlynniae_ASV1

trebouxia-pipeline --outdir results summarise \
  --ann-vcf results/06_variants/snpeff/variants.ann.vcf \
  --impact HIGH,MODERATE
```

---

## Functional annotation merge

```bash
python -m trebouxia_pipeline.merge_annotations \
  --gene-counts results/07_functional_summary/variant_counts_by_gene.tsv \
  --eggnog annotations/ASV1.emapper.annotations \
  --interpro annotations/ASV1.interproscan.tsv \
  --output results/07_functional_summary/variant_genes_functional_annotation.tsv
```

---

## Quick plots

```bash
python -m trebouxia_pipeline.plots \
  --input results/07_functional_summary/variant_effects_long.tsv \
  --column effect \
  --output results/07_functional_summary/snpEff_effects_barplot.png
```

---

## Recommended strict genomics additions before publication

1. Perform read-backed variant calling rather than relying only on assembly-to-assembly SNPs.
2. Mask repeats and low-complexity regions before interpreting private/exclusive genes.
3. Report callable genome size and coverage distribution for both strains.
4. Separate nuclear, chloroplast and mitochondrial comparisons.
5. For “exclusive genes”, require read coverage support and reciprocal protein searches.
6. Report whether non-synonymous variants are homozygous/fixed and avoid overinterpreting “high impact” calls without manual inspection.
7. Link key genes to phenotypes only as candidates unless validated by RNA-seq, proteomics or allele-specific assays.

---

## Associated Manuscript

This repository was developed as part of the research project described in the following manuscript:

> Blázquez M., Pizarro D., García-Muñoz A., Pino-Bodas R., de los Ríos A., Arzac M.I., García-Plazaola J.I., Fernández-Marín B., Artetxe U., Toledo-Gil R., Vallejo F., Casano L.M., Pérez-Vargas I., Pérez-Ortega S., & Gasulla F. Hidden strategies beneath a single nucleotide: contrasting physiology in two Trebouxia lynniae strains. Plant, Cell & Environment (currently under review).

The workflows implemented in this repository were developed to support the comparative genomic analyses performed in this study, including genome assembly assessment, whole-genome comparison, variant discovery, functional annotation, and interpretation of genetic differences between closely related Trebouxia strains. Please note that the manuscript is currently under peer review and some analyses, figures, and interpretations may evolve before publication.

## Project Status

Research Repository – Active Development. This repository accompanies an ongoing research project and is actively maintained. While the core analyses are fully reproducible, additional modules, documentation, and validation steps may be incorporated as the associated manuscript progresses through peer review. Users are welcome to explore, reproduce, and adapt the workflows for comparative genomics studies in microalgae, fungi, lichens, and other non-model organisms. CITE: https://doi.org/10.5281/zenodo.20700794

**Repository Maintainer:** David Pizarro, PhD  
Department of Pharmacology, Pharmacognosy and Botany  
Complutense University of Madrid (UCM)
