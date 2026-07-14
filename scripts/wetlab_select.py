"""
Selection stage of the wet-lab generation -> scoring -> selection pipeline.

Mirrors the old pipeline's "best" sets with a composite ranking per cell:
  1. keep valid, broad-AMP-positive candidates;
  2. for analogs, keep those in a "related but novel" similarity band to the
     parent (exclude near-duplicates and unrelated sequences);
  3. rank by predicted potency (APEX min MIC over the pathogen panel, ascending)
     then classifier consensus (descending);
  4. take the top K per cell.

APEX is expensive, so it is typically run only on a classifier shortlist. Where
apex.tsv is missing, this script still produces the classifier-ranked shortlist
and writes shortlist.fasta (the cheap input for score_wetlab_apex.sh); rerun
this script afterwards to get the APEX-ranked final selection.

NOTE: the new classifier ensemble has no hemolytic head (the old one did), so
the safety filter is omitted -- flagged in the output.

Usage:
  python scripts/wetlab_select.py [--wetlab data/wetlab] [--topk 100]
                                  [--pre 300] [--sim-min 40] [--sim-max 99]
"""
import argparse
from difflib import SequenceMatcher
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from wetlab_common import iter_cells, load_scored_cell


def similarity_to_parent(seq, parent):
    """Normalized sequence-similarity ratio (0-100) to the parent peptide."""
    return 100.0 * SequenceMatcher(None, seq, parent).ratio()


def select_cell(meta, topk, pre, sim_min, sim_max):
    df = load_scored_cell(meta["dir"])
    df = df[df["valid"]].copy()
    df["broad"] = pd.to_numeric(df["broad"], errors="coerce")
    df["consensus"] = pd.to_numeric(df["consensus"], errors="coerce")
    df["pathogen_min_mic"] = pd.to_numeric(df["pathogen_min_mic"], errors="coerce")

    # 1. AMP-positive candidates only
    cand = df[df["broad"] == 1].copy()
    if cand.empty:
        return None

    # 2. analog similarity band (related but novel)
    if meta["parent_seq"]:
        cand["similarity_to_parent"] = cand["Sequence"].apply(
            lambda s: similarity_to_parent(s, meta["parent_seq"]))
        cand = cand[(cand["similarity_to_parent"] >= sim_min) &
                    (cand["similarity_to_parent"] <= sim_max)]
    else:
        cand["similarity_to_parent"] = np.nan
    if cand.empty:
        return None

    cand = cand.drop_duplicates("Sequence")
    has_apex = cand["pathogen_min_mic"].notna().any()

    if has_apex:
        stage = "final"
        ranked = cand.sort_values(["pathogen_min_mic", "consensus"],
                                  ascending=[True, False]).head(topk)
    else:
        # classifier-only shortlist (needs APEX); broadest consensus first,
        # then closest-to-mid-band similarity for analogs
        stage = "needs_apex"
        cand["_simdist"] = (cand["similarity_to_parent"] - (sim_min + sim_max) / 2).abs()
        ranked = cand.sort_values(["consensus", "_simdist"],
                                  ascending=[False, True]).head(pre)
        ranked = ranked.drop(columns="_simdist")

    ranked.insert(0, "rank", range(1, len(ranked) + 1))
    for k in ("campaign", "parent", "mode", "sigma"):
        ranked[k] = meta[k]
    ranked["stage"] = stage
    ranked["hemolytic_filtered"] = False  # no hemolytic head in new ensemble
    return ranked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wetlab", type=Path, default=Path("data/wetlab"))
    ap.add_argument("--outdir", type=Path, default=Path("data/wetlab/selected"))
    ap.add_argument("--topk", type=int, default=100, help="final picks per scored cell")
    ap.add_argument("--pre", type=int, default=300, help="classifier shortlist size per unscored cell")
    ap.add_argument("--sim-min", type=float, default=40.0)
    ap.add_argument("--sim-max", type=float, default=99.0)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    keep = ["rank", "campaign", "parent", "mode", "sigma", "Sequence",
            "consensus", "broad", "pathogen_min_mic", "similarity_to_parent",
            "stage", "hemolytic_filtered"]
    all_rows, n_final, n_pre = [], 0, 0

    for meta in iter_cells(args.wetlab):
        sel = select_cell(meta, args.topk, args.pre, args.sim_min, args.sim_max)
        if sel is None or sel.empty:
            print(f"  {meta['cell']}: no candidates passed filters")
            continue
        sel = sel[[c for c in keep if c in sel.columns]]
        tag = meta["cell"].replace("/", "__")
        sel.to_csv(args.outdir / f"{tag}.csv", index=False)
        # shortlist FASTA for the cheap APEX run on cells still needing it
        if sel["stage"].iloc[0] == "needs_apex":
            n_pre += 1
            with open(meta["dir"] / "shortlist.fasta", "w") as fh:
                for _, r in sel.iterrows():
                    fh.write(f">{tag}_rank{r['rank']}\n{r['Sequence']}\n")
        else:
            n_final += 1
        all_rows.append(sel)

    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(args.wetlab / "selected_all.csv", index=False)
    print(f"\nWrote {args.wetlab/'selected_all.csv'}  ({len(combined)} rows)")
    print(f"  cells with APEX final selection: {n_final}")
    print(f"  cells with classifier shortlist (need APEX): {n_pre}  -> shortlist.fasta written")
    print("\n=== final-selection cells: top candidate per cell ===")
    fin = combined[combined["stage"] == "final"]
    if not fin.empty:
        best = fin[fin["rank"] == 1][["campaign", "parent", "mode", "sigma",
                                      "Sequence", "pathogen_min_mic", "consensus"]]
        print(best.to_string(index=False))


if __name__ == "__main__":
    main()
