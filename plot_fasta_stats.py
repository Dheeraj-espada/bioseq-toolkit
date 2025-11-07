#!/usr/bin/env python3
"""
plot_fasta_stats.py
Reads the CSV output from enhanced_seq_analyzer_cli.py and creates:
 - GC% distribution histogram (saved as gc_distribution.png)
 - Protein length distribution histogram (saved as protlen_distribution.png)

Usage:
    python3 plot_fasta_stats.py results.csv --out-dir plots
"""
import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

def plot_gc(df, out_path):
    plt.figure()
    df["gc_percent"].dropna().astype(float).plot.hist(bins=20)
    plt.xlabel("GC%")
    plt.ylabel("Count")
    plt.title("GC% distribution")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()

def plot_prot_len(df, out_path):
    plt.figure()
    # protein lengths may be string if empty; coerce to numeric, drop zeros/NaN
    prot = pd.to_numeric(df.get("prot_len", pd.Series()), errors="coerce").dropna()
    if len(prot) == 0:
        # fallback: try prot_seq_first50 length
        prot = pd.to_numeric(df.get("prot_len", pd.Series()), errors="coerce").dropna()
    plt.hist(prot, bins=20)
    plt.xlabel("Protein length (aa)")
    plt.ylabel("Count")
    plt.title("Protein length distribution")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()

def main():
    p = argparse.ArgumentParser(description="Plot FASTA/ORF stats from analyzer CSV")
    p.add_argument("csv", type=Path, help="CSV produced by enhanced_seq_analyzer_cli.py")
    p.add_argument("--out-dir", type=Path, default=Path("plots"), help="Output directory for plots")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.csv)

    gc_out = args.out_dir / "gc_distribution.png"
    prot_out = args.out_dir / "protlen_distribution.png"

    plot_gc(df, gc_out)
    plot_prot_len(df, prot_out)

    print("Wrote:", gc_out, prot_out)

if __name__ == "__main__":
    main()

