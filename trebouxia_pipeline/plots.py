#!/usr/bin/env python3
"""Simple plotting helpers for Trebouxia variant summaries."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a barplot from a TSV/CSV column")
    parser.add_argument("--input", required=True)
    parser.add_argument("--column", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top", type=int, default=30)
    args = parser.parse_args()

    sep = "\t" if str(args.input).endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(args.input, sep=sep)
    if args.column not in df.columns:
        raise ValueError(f"Column {args.column!r} not found. Available columns: {list(df.columns)}")

    counts = df[args.column].fillna("NA").value_counts().head(args.top)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    fig_height = max(4, 0.35 * len(counts))
    plt.figure(figsize=(8, fig_height))
    counts.sort_values().plot(kind="barh")
    plt.xlabel("Count")
    plt.ylabel(args.column)
    plt.title(f"Distribution of {args.column}")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
