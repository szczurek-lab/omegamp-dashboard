"""
Quality control for the wet-lab generated candidates (data/wetlab/).

One row per cell: generation sanity (counts, duplicates, validity, identity to
parent), scoring completeness/alignment, and headline score summaries
(AMP-positive rate, pathogen potency). Writes data/wetlab/qc_report.csv and
prints a per-campaign summary plus any flagged cells.

Usage: python scripts/wetlab_qc.py [--wetlab data/wetlab] [--out data/wetlab/qc_report.csv]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from wetlab_common import iter_cells, read_samples, load_scored_cell, scoring_complete


def qc_cell(meta):
    cell = meta["dir"]
    samples = read_samples(cell)
    seqs = [s for _, s in samples]
    n = len(seqs)
    uniq = len(set(seqs))
    lens = np.array([len(s) for s in seqs]) if seqs else np.array([0])
    row = {
        "campaign": meta["campaign"], "parent": meta["parent"],
        "mode": meta["mode"], "sigma": meta["sigma"], "cell": meta["cell"],
        "n": n, "n_unique": uniq,
        "dup_rate": round(1 - uniq / n, 4) if n else np.nan,
        "len_min": int(lens.min()), "len_med": int(np.median(lens)), "len_max": int(lens.max()),
        "scored": scoring_complete(cell),
    }
    # identity to parent: exact-sequence copies of the parent (analogs only)
    if meta["parent_seq"]:
        row["identical_to_parent_rate"] = round(np.mean([s == meta["parent_seq"] for s in seqs]), 4)
    else:
        row["identical_to_parent_rate"] = np.nan

    if row["scored"]:
        df = load_scored_cell(cell)
        row["rows_aligned"] = bool(df["pathogen_min_mic"].notna().any() and len(df) == n)
        row["valid_rate"] = round(df["valid"].mean(), 4)
        row["broad_amp_rate"] = round(pd.to_numeric(df["broad"], errors="coerce").mean(), 4)
        mic = pd.to_numeric(df["pathogen_min_mic"], errors="coerce")
        row["median_pathogen_mic"] = round(float(mic.median()), 2)
        row["pct_mic_le8"] = round(float((mic <= 8).mean()) * 100, 1)
        row["mic_missing"] = int(mic.isna().sum())
    else:
        valid_rate = np.mean([len(s) > 0 and set(s) <= set("ACDEFGHIKLMNPQRSTVWY") for s in seqs]) if seqs else 0
        row.update({"rows_aligned": False, "valid_rate": round(valid_rate, 4),
                    "broad_amp_rate": np.nan, "median_pathogen_mic": np.nan,
                    "pct_mic_le8": np.nan, "mic_missing": np.nan})
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wetlab", type=Path, default=Path("data/wetlab"))
    ap.add_argument("--out", type=Path, default=Path("data/wetlab/qc_report.csv"))
    args = ap.parse_args()

    rows = [qc_cell(m) for m in iter_cells(args.wetlab)]
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)

    print(f"Wrote {args.out}  ({len(df)} cells)\n")
    print("=== per-campaign summary ===")
    agg = df.groupby("campaign").agg(
        cells=("cell", "size"), scored=("scored", "sum"),
        total_seqs=("n", "sum"), mean_dup=("dup_rate", "mean"),
        mean_valid=("valid_rate", "mean"), mean_amp_rate=("broad_amp_rate", "mean"),
        med_pathogen_mic=("median_pathogen_mic", "median"),
    ).round(3)
    print(agg.to_string())

    print("\n=== flags ===")
    flagged = False
    for _, r in df.iterrows():
        issues = []
        if not r["scored"]:
            issues.append("UNSCORED")
        elif not r["rows_aligned"]:
            issues.append("ROWS_NOT_ALIGNED")
        if pd.notna(r["valid_rate"]) and r["valid_rate"] < 0.99:
            issues.append(f"valid={r['valid_rate']}")
        if pd.notna(r["dup_rate"]) and r["dup_rate"] > 0.5:
            issues.append(f"dup={r['dup_rate']}")
        if pd.notna(r["identical_to_parent_rate"]) and r["identical_to_parent_rate"] > 0.1:
            issues.append(f"parent_copies={r['identical_to_parent_rate']}")
        if pd.notna(r.get("mic_missing")) and r["mic_missing"] > 0:
            issues.append(f"mic_missing={int(r['mic_missing'])}")
        if issues:
            flagged = True
            print(f"  {r['cell']}: {', '.join(issues)}")
    if not flagged:
        print("  none")


if __name__ == "__main__":
    main()
