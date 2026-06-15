#!/usr/bin/env python3
"""Merge variant-impact summaries with EggNOG-mapper and InterProScan annotations."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def read_eggnog(path: str | Path) -> pd.DataFrame:
    rows = []
    header = None
    with open(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header and len(fields) == len(header):
                rows.append(dict(zip(header, fields)))
            else:
                # EggNOG classic columns fallback.
                rows.append({"query": fields[0], "eggNOG_raw": "\t".join(fields[1:])})
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["gene_id"])
    if "query" in df.columns:
        df = df.rename(columns={"query": "gene_id"})
    elif "#query" in df.columns:
        df = df.rename(columns={"#query": "gene_id"})
    keep = [c for c in ["gene_id", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "PFAMs", "eggNOG_raw"] if c in df.columns]
    return df[keep].drop_duplicates("gene_id")


def read_interpro(path: str | Path) -> pd.DataFrame:
    # Standard InterProScan TSV columns, no header.
    cols = [
        "gene_id", "md5", "length", "analysis", "signature_accession",
        "signature_description", "start", "end", "score", "status", "date",
        "interpro_accession", "interpro_description", "go_terms", "pathways",
    ]
    df = pd.read_csv(path, sep="\t", header=None, comment="#", names=cols, usecols=range(15))
    if df.empty:
        return pd.DataFrame(columns=["gene_id"])
    grouped = df.groupby("gene_id", as_index=False).agg({
        "analysis": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
        "signature_accession": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
        "signature_description": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
        "interpro_accession": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
        "interpro_description": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
        "go_terms": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
        "pathways": lambda x: ";".join(sorted(set(map(str, x.dropna())))),
    })
    return grouped


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge variant counts with functional annotations")
    parser.add_argument("--gene-counts", required=True)
    parser.add_argument("--eggnog", required=True)
    parser.add_argument("--interpro", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    counts = pd.read_csv(args.gene_counts, sep="\t")
    eggnog = read_eggnog(args.eggnog)
    interpro = read_interpro(args.interpro)

    merged = counts.merge(eggnog, on="gene_id", how="left")
    merged = merged.merge(interpro, on="gene_id", how="left", suffixes=("_eggnog", "_interpro"))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
