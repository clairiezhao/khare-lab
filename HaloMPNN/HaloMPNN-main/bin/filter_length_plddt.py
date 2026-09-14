#!/usr/bin/env python3
"""
filter_length_plddt.py  <min_length> <min_plddt>
                        --input  metadata.csv
                        --output metadata_filtered.csv
                        --stats  filter_stats.txt

Filters a metadata CSV by sequence length and mean pLDDT, saves stats,
and produces histograms for both columns (before and after filtering).
"""

import argparse
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("min_length", type=int,
                   help="Minimum sequence length (inclusive)")
    p.add_argument("min_plddt", type=float,
                   help="Minimum mean pLDDT (inclusive)")
    p.add_argument("--input",  required=True, help="Input CSV file")
    p.add_argument("--output", required=True, help="Filtered output CSV file")
    p.add_argument("--stats",  required=True, help="Plain-text stats file")
    return p.parse_args()


# ── Plotting helper ────────────────────────────────────────────────────────────
def plot_hist(series, title, xlabel, ylabel, outpath, bins=20, vline=None):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(series.dropna(), bins=bins, edgecolor="white", linewidth=0.4)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    if vline is not None:
        ax.axvline(vline, color="red", linestyle="--", linewidth=1.2,
                   label=f"threshold = {vline}")
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()

    df_raw = pd.read_csv(args.input)
    n_start = len(df_raw)

    # ── Step 1: filter by length ───────────────────────────────────────────────
    mask_len = df_raw["length"] >= args.min_length
    n_fail_len = (~mask_len).sum()
    df_after_len = df_raw[mask_len].copy()

    # ── Step 2: filter by pLDDT ───────────────────────────────────────────────
    mask_plddt = df_after_len["mean_plddt"] >= args.min_plddt
    n_fail_plddt = (~mask_plddt).sum()
    df_filtered = df_after_len[mask_plddt].copy()

    n_final = len(df_filtered)
    n_total_removed = n_start - n_final

    # ── Save filtered CSV ─────────────────────────────────────────────────────
    df_filtered.to_csv(args.output, index=False)

    # ── Save stats ────────────────────────────────────────────────────────────
    stats_lines = [
        f"Input sequences:                  {n_start}",
        f"Removed (length < {args.min_length}):           {n_fail_len}",
        f"Remaining after length filter:    {n_start - n_fail_len}",
        f"Removed (mean_plddt < {args.min_plddt}):      {n_fail_plddt}",
        f"Remaining after pLDDT filter:     {n_final}",
        f"Total removed:                    {n_total_removed}",
        f"Total kept:                       {n_final}",
    ]
    stats_text = "\n".join(stats_lines)
    print(stats_text)
    with open(args.stats, "w") as fh:
        fh.write(stats_text + "\n")

    # ── Histograms (before filtering) ─────────────────────────────────────────
    plot_hist(
        series   = df_raw["length"],
        title    = f"Sequence length — before filtering  (n={n_start})",
        xlabel   = "Sequence length (aa)",
        ylabel   = "Frequency",
        outpath  = "hist_length_before.png",
        vline    = args.min_length,
    )
    plot_hist(
        series   = df_raw["mean_plddt"],
        title    = f"Mean pLDDT — before filtering  (n={n_start})",
        xlabel   = "Mean pLDDT",
        ylabel   = "Frequency",
        outpath  = "hist_plddt_before.png",
        vline    = args.min_plddt,
    )

    # ── Histograms (after filtering) ──────────────────────────────────────────
    plot_hist(
        series   = df_filtered["length"],
        title    = f"Sequence length — after filtering  (n={n_final})",
        xlabel   = "Sequence length (aa)",
        ylabel   = "Frequency",
        outpath  = "hist_length_after.png",
    )
    plot_hist(
        series   = df_filtered["mean_plddt"],
        title    = f"Mean pLDDT — after filtering  (n={n_final})",
        xlabel   = "Mean pLDDT",
        ylabel   = "Frequency",
        outpath  = "hist_plddt_after.png",
    )

    print(f"\nSaved: {args.output}, {args.stats}")
    print("Saved: hist_length_before.png, hist_plddt_before.png")
    print("Saved: hist_length_after.png,  hist_plddt_after.png")


if __name__ == "__main__":
    main()
