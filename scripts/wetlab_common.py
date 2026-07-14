"""
Shared helpers for the wet-lab generation -> scoring -> selection pipeline.

A wet-lab "cell" is one leaf directory under ``data/wetlab/`` holding:
  samples.fasta     generated candidates (>sequence_N)
  predictions.csv   classifier ensemble, BINARY 0/1 columns
                    (broad-classifier + 5 species + 7 strains)
  apex.tsv          APEX MIC (uM) per strain, 34 columns, indexed by sequence

Campaigns: denovo/ (modes), analog/<parent>/<mode>[/sigma*], motif/<motif>[/sigma*].

Potency is scored against a curated PATHOGEN panel only -- the APEX panel also
contains gut commensals (Bacteroides, Akkermansia, ...) and the probiotic
E. coli Nissle, where killing is NOT desirable, so those are excluded.
"""
from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO

# 15 clinically relevant pathogens (ESKAPE + MRSA/VRE + Salmonella + Listeria).
# Everything else in the 34-strain APEX panel is a gut commensal or the
# probiotic E. coli Nissle and is excluded from potency scoring.
PATHOGEN_STRAINS = [
    "E. coli ATCC11775", "P. aeruginosa PAO1", "P. aeruginosa PA14",
    "S. aureus ATCC12600", "E. coli AIG221", "E. coli AIG222",
    "K. pneumoniae ATCC13883", "A. baumannii ATCC19606",
    "S. aureus (ATCC BAA-1556) - MRSA",
    "vancomycin-resistant E. faecalis ATCC700802",
    "vancomycin-resistant E. faecium ATCC700221",
    "Salmonella enterica ATCC 9150 (BEIRES NR-515)",
    "Salmonella enterica (BEIRES NR-170)",
    "Salmonella enterica ATCC 9150 (BEIRES NR-174)",
    "L. monocytogenes ATCC 19111 (BEIRES NR-106)",
]

CLASSIFIER_COLS = None  # resolved per-file (all *-classifier columns)
AA = set("ACDEFGHIKLMNPQRSTVWY")


def iter_cells(wetlab_root):
    """Yield dicts describing every cell (dir with samples.fasta) under root."""
    root = Path(wetlab_root)
    for fasta in sorted(root.rglob("samples.fasta")):
        cell = fasta.parent
        rel = cell.relative_to(root)
        parts = rel.parts
        campaign = parts[0]
        meta = {"campaign": campaign, "cell": str(rel), "dir": cell,
                "parent": None, "mode": None, "sigma": None, "parent_seq": None}
        if campaign == "analog":
            meta["parent"] = parts[1]
            meta["mode"] = parts[2] if len(parts) > 2 else None
            meta["parent_seq"] = _read_parent(cell)
        elif campaign == "motif":
            meta["parent"] = parts[1]            # motif name
            meta["mode"] = "motif"
        elif campaign == "denovo":
            meta["mode"] = parts[1] if len(parts) > 1 else None
        sig = re.search(r"sigma([\d.]+)", str(rel))
        if sig:
            meta["sigma"] = float(sig.group(1))
        yield meta


def _read_parent(cell):
    """Walk up from a cell dir to find parent.fasta; return its sequence."""
    for d in [cell, *cell.parents]:
        pf = d / "parent.fasta"
        if pf.exists():
            recs = list(SeqIO.parse(str(pf), "fasta"))
            if recs:
                return str(recs[0].seq)
        if d.name == "analog":
            break
    return None


def read_samples(cell_dir):
    """Return list of (id, seq) from samples.fasta in FASTA order."""
    return [(r.id, str(r.seq)) for r in SeqIO.parse(str(Path(cell_dir) / "samples.fasta"), "fasta")]


def load_scored_cell(cell_dir):
    """Merge classifier + APEX onto each generated sequence.

    Returns a DataFrame with one row per generated sequence:
      Sequence, valid, consensus (# positive classifiers), broad (0/1),
      pathogen_min_mic (min APEX uM over PATHOGEN_STRAINS), and the raw
      classifier columns. Missing score files -> those columns are NaN.
    """
    cell = Path(cell_dir)
    samples = read_samples(cell)
    df = pd.DataFrame(samples, columns=["Id", "Sequence"])
    df["valid"] = df["Sequence"].apply(lambda s: len(s) > 0 and set(s) <= AA)

    pred = cell / "predictions.csv"
    if pred.exists():
        cdf = pd.read_csv(pred)
        clf_cols = [c for c in cdf.columns if c.endswith("-classifier")]
        cons = cdf.set_index("Sequence")[clf_cols].groupby(level=0).max()
        df["consensus"] = df["Sequence"].map(cons.sum(axis=1))
        if "broad-classifier" in clf_cols:
            df["broad"] = df["Sequence"].map(cons["broad-classifier"])
    else:
        df["consensus"] = pd.NA
        df["broad"] = pd.NA

    apex = cell / "apex.tsv"
    if apex.exists():
        adf = pd.read_csv(apex, sep="\t", index_col=0)
        cols = [c for c in PATHOGEN_STRAINS if c in adf.columns]
        pmic = adf[cols].min(axis=1).groupby(level=0).min()
        df["pathogen_min_mic"] = df["Sequence"].map(pmic)
    else:
        df["pathogen_min_mic"] = pd.NA

    return df


def scoring_complete(cell_dir):
    cell = Path(cell_dir)
    return (cell / "predictions.csv").exists() and (cell / "apex.tsv").exists()
