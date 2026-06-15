#!/usr/bin/env python3
"""Command-line interface for the Trebouxia comparative genomics pipeline.

The CLI creates transparent shell commands for common comparative-genomics
steps. By default it executes them. Use --dry-run to print commands without
running them.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shlex
import subprocess
from pathlib import Path
from typing import Iterable, List, Dict


def run_command(cmd: str, dry_run: bool = False) -> None:
    print(f"[trebouxia-pipeline] {cmd}")
    if not dry_run:
        subprocess.run(cmd, shell=True, check=True)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_samples(samples_csv: str | Path) -> List[Dict[str, str]]:
    with open(samples_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sample", "r1", "r2"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in samples.csv: {missing}")
        return list(reader)


def sample_by_name(samples: List[Dict[str, str]], name: str) -> Dict[str, str]:
    for sample in samples:
        if sample["sample"] == name:
            return sample
    raise ValueError(f"Sample {name!r} not found in samplesheet")


def cmd_qc(args: argparse.Namespace) -> None:
    samples = read_samples(args.samples)
    raw_qc = Path(args.outdir) / "01_qc" / "raw_fastqc"
    multiqc = Path(args.outdir) / "01_qc" / "multiqc_raw"
    ensure_dir(raw_qc)
    ensure_dir(multiqc)

    fastqs = []
    for s in samples:
        fastqs.extend([s["r1"], s["r2"]])
    run_command(f"fastqc -t {args.threads} -o {raw_qc} " + " ".join(map(shlex.quote, fastqs)), args.dry_run)
    run_command(f"multiqc -o {multiqc} {raw_qc}", args.dry_run)


def cmd_trim(args: argparse.Namespace) -> None:
    samples = read_samples(args.samples)
    out = Path(args.outdir) / "01_qc" / "trimmed"
    qc = Path(args.outdir) / "01_qc" / "fastp_reports"
    ensure_dir(out)
    ensure_dir(qc)

    for s in samples:
        sample = s["sample"]
        r1_out = out / f"{sample}_R1.trimmed.fastq.gz"
        r2_out = out / f"{sample}_R2.trimmed.fastq.gz"
        html = qc / f"{sample}.fastp.html"
        json = qc / f"{sample}.fastp.json"
        cmd = (
            f"fastp -w {args.threads} "
            f"-i {shlex.quote(s['r1'])} -I {shlex.quote(s['r2'])} "
            f"-o {r1_out} -O {r2_out} "
            f"--html {html} --json {json}"
        )
        run_command(cmd, args.dry_run)


def cmd_assemble(args: argparse.Namespace) -> None:
    samples = read_samples(args.samples)
    trimmed = Path(args.outdir) / "01_qc" / "trimmed"
    outroot = Path(args.outdir) / "02_assembly"
    ensure_dir(outroot)

    for s in samples:
        sample = s["sample"]
        r1 = trimmed / f"{sample}_R1.trimmed.fastq.gz"
        r2 = trimmed / f"{sample}_R2.trimmed.fastq.gz"
        if not r1.exists():
            r1 = Path(s["r1"])
            r2 = Path(s["r2"])
        out = outroot / sample
        ensure_dir(out)
        cmd = f"spades.py -1 {r1} -2 {r2} -o {out} -t {args.threads} --careful"
        run_command(cmd, args.dry_run)


def cmd_assembly_qc(args: argparse.Namespace) -> None:
    samples = read_samples(args.samples)
    outroot = Path(args.outdir)
    quast_dir = outroot / "03_assembly_qc" / "quast"
    busco_dir = outroot / "03_assembly_qc" / "busco"
    ensure_dir(quast_dir)
    ensure_dir(busco_dir)

    assemblies = []
    for s in samples:
        asm = outroot / "02_assembly" / s["sample"] / "scaffolds.fasta"
        assemblies.append(str(asm))
        run_command(
            f"busco -i {asm} -o {s['sample']} -m genome -l {args.busco_lineage} -c {args.threads} --out_path {busco_dir}",
            args.dry_run,
        )
    run_command(f"quast.py -t {args.threads} -o {quast_dir} " + " ".join(assemblies), args.dry_run)


def cmd_compare(args: argparse.Namespace) -> None:
    samples = read_samples(args.samples)
    if len(samples) != 2:
        raise ValueError("The compare command expects exactly two samples")
    outroot = Path(args.outdir)
    a = samples[0]["sample"]
    b = samples[1]["sample"]
    asm_a = outroot / "02_assembly" / a / "scaffolds.fasta"
    asm_b = outroot / "02_assembly" / b / "scaffolds.fasta"
    out = outroot / "04_genome_comparison"
    ensure_dir(out)

    run_command(f"fastANI -q {asm_a} -r {asm_b} -o {out}/{a}_vs_{b}.fastani.tsv -t {args.threads}", args.dry_run)
    run_command(f"nucmer --threads={args.threads} --prefix={out}/{a}_vs_{b} {asm_a} {asm_b}", args.dry_run)
    run_command(f"delta-filter -1 {out}/{a}_vs_{b}.delta > {out}/{a}_vs_{b}.1delta", args.dry_run)
    run_command(f"show-coords -rcl {out}/{a}_vs_{b}.1delta > {out}/{a}_vs_{b}.coords", args.dry_run)
    run_command(f"show-snps -Clr {out}/{a}_vs_{b}.1delta > {out}/{a}_vs_{b}.snps", args.dry_run)


def cmd_variants(args: argparse.Namespace) -> None:
    samples = read_samples(args.samples)
    query = sample_by_name(samples, args.query_sample)
    outroot = Path(args.outdir)
    out = outroot / "06_variants"
    bamdir = out / "bam"
    ensure_dir(out)
    ensure_dir(bamdir)
    reference = Path(args.reference)
    sample = args.query_sample

    r1 = outroot / "01_qc" / "trimmed" / f"{sample}_R1.trimmed.fastq.gz"
    r2 = outroot / "01_qc" / "trimmed" / f"{sample}_R2.trimmed.fastq.gz"
    if not r1.exists():
        r1 = Path(query["r1"])
        r2 = Path(query["r2"])

    bam = bamdir / f"{sample}_vs_reference.sorted.bam"
    raw_vcf = out / f"{sample}_vs_ASV1.raw.vcf.gz"
    filt_vcf = out / f"{sample}_vs_ASV1.filtered.vcf.gz"

    run_command(f"bwa index {reference}", args.dry_run)
    run_command(
        f"bwa mem -t {args.threads} {reference} {r1} {r2} | samtools sort -@ {args.threads} -o {bam}",
        args.dry_run,
    )
    run_command(f"samtools index {bam}", args.dry_run)
    run_command(f"samtools coverage {bam} > {out}/{sample}_coverage.tsv", args.dry_run)
    run_command(
        f"bcftools mpileup -Ou -f {reference} {bam} | bcftools call -mv -Oz -o {raw_vcf}",
        args.dry_run,
    )
    run_command(f"bcftools index {raw_vcf}", args.dry_run)
    run_command(
        f"bcftools filter -Oz -o {filt_vcf} -i 'QUAL>={args.min_qual} && DP>={args.min_depth}' {raw_vcf}",
        args.dry_run,
    )
    run_command(f"bcftools index {filt_vcf}", args.dry_run)


def cmd_snpeff(args: argparse.Namespace) -> None:
    out = Path(args.outdir) / "06_variants" / "snpeff"
    ensure_dir(out)
    ann = out / "variants.ann.vcf"
    html = out / "snpEff_summary.html"
    cmd = f"snpEff -stats {html} {shlex.quote(args.snpeff_db)} {shlex.quote(args.vcf)} > {ann}"
    run_command(cmd, args.dry_run)


def parse_ann_field(info: str) -> List[Dict[str, str]]:
    """Parse SnpEff ANN records from the INFO column."""
    ann_values = []
    for part in info.split(";"):
        if part.startswith("ANN="):
            ann_values = part.replace("ANN=", "", 1).split(",")
            break
    parsed = []
    keys = [
        "allele", "effect", "impact", "gene_name", "gene_id", "feature_type",
        "feature_id", "transcript_biotype", "rank", "hgvs_c", "hgvs_p",
        "cdna_pos", "cds_pos", "aa_pos", "distance", "errors",
    ]
    for value in ann_values:
        fields = value.split("|")
        fields += [""] * (len(keys) - len(fields))
        parsed.append(dict(zip(keys, fields)))
    return parsed


def cmd_summarise(args: argparse.Namespace) -> None:
    out = Path(args.outdir) / "07_functional_summary"
    ensure_dir(out)
    impacts = {x.strip() for x in args.impact.split(",") if x.strip()}
    long_rows = []
    counts = {}

    opener = gzip.open if str(args.ann_vcf).endswith(".gz") else open
    with opener(args.ann_vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = fields[:8]
            for ann in parse_ann_field(info):
                if impacts and ann["impact"] not in impacts:
                    continue
                gene = ann.get("gene_id") or ann.get("gene_name") or "NA"
                row = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "filter": flt,
                    **ann,
                }
                long_rows.append(row)
                counts.setdefault(gene, {"gene_id": gene, "n_variants": 0, "HIGH": 0, "MODERATE": 0, "LOW": 0, "MODIFIER": 0})
                counts[gene]["n_variants"] += 1
                if ann["impact"] in counts[gene]:
                    counts[gene][ann["impact"]] += 1

    import pandas as pd

    pd.DataFrame(long_rows).to_csv(out / "variant_effects_long.tsv", sep="\t", index=False)
    pd.DataFrame(counts.values()).sort_values("n_variants", ascending=False).to_csv(
        out / "variant_counts_by_gene.tsv", sep="\t", index=False
    )
    print(f"Wrote {out / 'variant_effects_long.tsv'}")
    print(f"Wrote {out / 'variant_counts_by_gene.tsv'}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Trebouxia comparative genomics CLI")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them")

    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--samples", required=True, help="CSV file with sample,r1,r2 columns")
    common.add_argument("--threads", type=int, default=8)

    p = sub.add_parser("qc", parents=[common], help="Run FastQC and MultiQC on raw FASTQ files")
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("trim", parents=[common], help="Trim reads with fastp")
    p.set_defaults(func=cmd_trim)

    p = sub.add_parser("assemble", parents=[common], help="Assemble genomes with SPAdes")
    p.set_defaults(func=cmd_assemble)

    p = sub.add_parser("assembly-qc", parents=[common], help="Run QUAST and BUSCO")
    p.add_argument("--busco-lineage", default="chlorophyta_odb10")
    p.set_defaults(func=cmd_assembly_qc)

    p = sub.add_parser("compare", parents=[common], help="Run FastANI and MUMmer comparison")
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("variants", parents=[common], help="Map reads and call variants with bcftools")
    p.add_argument("--reference", required=True, help="Reference FASTA, e.g. ASV1 assembly")
    p.add_argument("--query-sample", required=True, help="Sample to map against the reference")
    p.add_argument("--min-qual", type=float, default=30)
    p.add_argument("--min-depth", type=int, default=10)
    p.set_defaults(func=cmd_variants)

    p = sub.add_parser("snpeff", help="Annotate VCF with SnpEff")
    p.add_argument("--vcf", required=True)
    p.add_argument("--snpeff-db", required=True)
    p.set_defaults(func=cmd_snpeff)

    p = sub.add_parser("summarise", help="Summarise SnpEff ANN effects by gene")
    p.add_argument("--ann-vcf", required=True)
    p.add_argument("--impact", default="HIGH,MODERATE")
    p.set_defaults(func=cmd_summarise)

    return parser


def main(argv: Iterable[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
