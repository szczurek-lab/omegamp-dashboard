#!/usr/bin/env python3
"""Build processed dashboard CSVs from Marcelo's raw data files.

Run from the repository root:
    python scripts/build_data.py

Override the raw-data directory:
    OMEGAMP_RAW=/path/to/raw python scripts/build_data.py

Inputs  (read-only): $OMEGAMP_RAW/  (default: data/raw/)
                     data/omegamp_reference_table.csv
Outputs (overwrite): data/*.csv

Censoring rules
---------------
- MIC:        NaN (inactive) → ">64"
- HC50/CC50:  value > 128 or non-finite → ">128"

Known ID collision
------------------
sequence_9310 appears in both the de novo and bZIP-analog source files with
different sequences.  The reference table disambiguates via the alias
sequence_9310b → Ω-MT-bZIP-4.  Files in a bZIP/analog context apply this
override automatically; all other contexts map sequence_9310 → Ω-DP-68.
"""

import csv as csv_module
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path(os.environ.get("OMEGAMP_RAW", "data/raw"))
OUT = Path("data")
REF_TABLE = OUT / "omegamp_reference_table.csv"

# ── Canonical MIC pathogen column order ──────────────────────────────────────
MIC_PATHOGENS = [
    "A. baumannii ATCC 19606 (-)",
    "A. baumannii ATCC BAA-1605 (-) - CGTPACCIMRA",
    "E. cloacae ATCC 13047 (-)",
    "E. coli ATCC 11775 (-)",
    "E. coli AIC221 (-)",
    "E. coli AIC222 - CRE (-)",
    "E. coli ATCC BAA-3170 (-) - CRE",
    "K. pneumoniae ATCC 13883 (-)",
    "K. pneumoniae ATCC BAA-2342 (-) - EIRK",
    "P. aeruginosa PAO1 (-)",
    "P. aeruginosa PA14 (-)",
    "P. aeruginosa ATCC BAA-3197 (-) - FBCRP",
    "S. enterica ATCC 9150 (-)",
    "S. enterica Typhimurtium ATCC 700720",
    "B. subtilis ATCC 23857 (+)",
    "S. aureus ATCC 12600 (+)",
    "S. aureus ATCC BAA-1556 - MRSA (+)",
    "E. faecalis ATCC 700802 - VRE (+)",
    "E. faecium ATCC 700221 - VRE (+)",
    "E. coli K-12 BW25113 (-)",
]

# ── LPS binding ───────────────────────────────────────────────────────────────
# family key (prefix of _at_45min*.txt filenames) → canonical family name
LPS_FAMILIES = {
    "Cecropin-P1": "cecropin",
    "LG21":        "LG21",
    "pa4":         "pa4",
    "sarcotoxin-1A": "sarcotoxin",
}
# Prototype/control name normalisation within LPS files
LPS_NAME_NORM = {
    "cecropin-P1":   "cecropin",
    "LG21":          "LG21",
    "pa4":           "pa4",
    "sarcotoxin-1A": "sarcotoxin",
}

# ── BeStSel secondary structure ───────────────────────────────────────────────
# Raw file groups columns as 5 fractions × 4 solvents (H2O, TFE/H2O, SDS/H2O, MeOH/H2O).
# The previous repo bestsel.csv had TFE and MeOH labels swapped; this script
# assigns names correctly from the raw file's annotation.
BESTSEL_FRACTIONS = ["fH", "f_beta_anti", "f_beta_par", "fturn", "fothers"]
BESTSEL_SOLVENTS  = ["H2O", "TFE_H2O", "SDS_H2O", "MeOH_H2O"]   # raw file order
BESTSEL_COLS = [f"{f}_{s}" for f in BESTSEL_FRACTIONS for s in BESTSEL_SOLVENTS]

# Every de novo peptide has measured BeStSel data. If bestsel.csv drops below
# this, the per-solvent raw CSVs are missing and figure S4 / the dashboard
# silently lose peptides — fail the build instead (see 25-vs-95 regression).
EXPECTED_DENOVO_BESTSEL = 95

# Per-solvent BeStSel CSV suffixes, in BESTSEL_SOLVENTS order.
_SS_SOLVENT_SUFFIXES = [
    "_BeStSel_H2O_heatmap.csv",
    "_BeStSel_TFE_H2O_3_2_v_v_heatmap.csv",
    "_BeStSel_SDS_10mM_H2O_heatmap.csv",
    "_BeStSel_MeOH_H2O_1_1_v_v_heatmap.csv",
]
_SS_PER_SOLVENT_RE = re.compile(
    r"_BeStSel_(?:H2O|TFE_H2O_3_2_v_v|SDS_10mM_H2O|MeOH_H2O_1_1_v_v)_heatmap\.csv$"
)

# ── sequence_9310 collision override (bZIP / analog context) ─────────────────
BZIP_OVERRIDE = {"sequence_9310": "sequence_9310b"}


# ─────────────────────────────────────────────────────────────────────────────
# Reference-table helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_ref():
    """Return (by_seqid, by_seq) — both map to short_name."""
    df = pd.read_csv(REF_TABLE)
    by_seqid = (
        df.dropna(subset=["sequence_id"])
        .set_index("sequence_id")["short_name"]
        .to_dict()
    )
    by_seq = (
        df.dropna(subset=["sequence"])
        .set_index("sequence")["short_name"]
        .to_dict()
    )
    return by_seqid, by_seq


def resolve(pid: str, seq: str, by_seqid: dict, by_seq: dict,
            overrides: dict | None = None) -> str:
    """Map raw peptide_id (+ optional amino-acid sequence) to short_name.

    Sequence lookup takes priority over ID lookup when a sequence is provided:
    the amino-acid sequence is unique, while a peptide_id may be ambiguous
    (e.g. sequence_9310 appears in both the de-novo and bZIP-analog source files
    with different sequences).

    overrides: local id→id substitutions applied before ID lookup (for cases
               where no sequence is available but the context is known).
    Falls back to pid as-is for controls / prototypes not in the reference table.
    """
    # 1. Sequence is the most specific key — use it first
    if seq and seq in by_seq:
        return by_seq[seq]
    # 2. Apply context-specific ID override before the general lookup
    effective = overrides.get(pid, pid) if overrides else pid
    if effective in by_seqid:
        return by_seqid[effective]
    # 3. Fallback: return pid as-is (controls, templates, benchmarks)
    return pid


# ─────────────────────────────────────────────────────────────────────────────
# Censoring helpers
# ─────────────────────────────────────────────────────────────────────────────

def censor_mic(val) -> str:
    """NaN → '>64'; numeric values pass through as strings."""
    if pd.isna(val):
        return ">64"
    # Keep integer-looking values without trailing .0
    f = float(val)
    return str(int(f)) if f == int(f) else str(f)


def censor_128(val):
    """Non-finite or >128 → '>128'; others pass through as float."""
    try:
        v = float(val)
    except (TypeError, ValueError):
        return val
    if not np.isfinite(v) or v > 128.0:
        return ">128"
    return v


# ─────────────────────────────────────────────────────────────────────────────
# Build functions
# ─────────────────────────────────────────────────────────────────────────────

def build_mic(by_seqid, by_seq):
    """
    Sources
    -------
    Summary_pathogens.csv           1st test batch (de novo); has Peptide+Sequence columns.
    Summary_pathogens_2nd_set.csv   2nd test batch (analogs + more); blank first column.
    MIC_Summary_Pathogens_with_lower_than_1uM_results.csv
                                    Re-reads, below 1 uM, of peptides that bottomed
                                    out at MIC = 1 uM in the batches above.  Its
                                    header row is misaligned with its data and is
                                    corrected on load -- see the note further down.

    In the 2nd-set file sequence_9310 is the bZIP-4 analog → apply BZIP_OVERRIDE.
    """
    # ── 1st set (Peptide + Sequence → use sequence for disambiguation) ────────
    df1 = pd.read_csv(RAW / "mic" / "Summary_pathogens.csv")
    df1["short_name"] = [
        resolve(p, str(s) if pd.notna(s) else "", by_seqid, by_seq)
        for p, s in zip(df1["Peptide"], df1["Sequence"])
    ]
    df1 = df1.drop(columns=["Peptide", "Sequence"])

    # ── 2nd set (sequence_id only) ────────────────────────────────────────────
    # Special case: "sequence_9310" appears TWICE in this file with no Sequence
    # column to distinguish the entries.  By position:
    #   1st occurrence → Ω-DP-68       (de-novo, low MIC values)
    #   2nd occurrence → Ω-MT-bZIP-4   (bZIP analog)
    df2 = pd.read_csv(RAW / "mic" / "Summary_pathogens_2nd_set.csv")
    df2 = df2.rename(columns={df2.columns[0]: "pid"})
    seen_9310 = 0
    short_names_2nd = []
    for pid in df2["pid"]:
        if pid == "sequence_9310":
            seen_9310 += 1
            short_names_2nd.append("Ω-DP-68" if seen_9310 == 1 else "Ω-MT-bZIP-4")
        else:
            short_names_2nd.append(resolve(pid, "", by_seqid, by_seq))
    df2["short_name"] = short_names_2nd
    df2 = df2.drop(columns=["pid"])

    # ── Merge: align to canonical column order, concatenate ───────────────────
    combined = pd.concat(
        [
            df1.reindex(columns=["short_name"] + MIC_PATHOGENS),
            df2.reindex(columns=["short_name"] + MIC_PATHOGENS),
        ],
        ignore_index=True,
    )

    # ── Deduplicate: for conflicting short_names keep the more complete row ───
    # Sort so the row with fewest NaNs comes first, then deduplicate keeping first.
    combined["_n_null"] = combined[MIC_PATHOGENS].isnull().sum(axis=1)
    combined = (
        combined
        .sort_values("_n_null")
        .drop_duplicates(subset=["short_name"], keep="first")
        .drop(columns=["_n_null"])
        .reset_index(drop=True)
    )

    # ── 3rd source: sub-1 uM re-reads ────────────────────────────────────────
    # Peptides that bottomed out at MIC = 1 uM were re-measured below 1 uM.
    #
    # The shipped file's HEADER ROW IS OUT OF STEP WITH ITS DATA.  Its data
    # columns follow the 1st-set layout (E. coli K-12 BW25113 sits 8th, right
    # after E. coli BAA-3170); the shipped header is the 2nd-set layout (K-12
    # last).  Read as shipped, every value from K. pneumoniae ATCC 13883 onward
    # lands on the strain one column to its left.  So we discard the shipped
    # header and relabel the data with the 1st-set layout.
    #
    # Verified two independent ways: under this layout all 25 first-set peptides
    # match Summary_pathogens.csv exactly, and every peptide with no sub-1 value
    # matches Summary_pathogens_2nd_set.csv exactly (99/99).  The resulting merge
    # touches only sub-1 cells -- no strain reassignment, no activity-threshold
    # crossings.
    #
    # The Polymyxin B / Levofloxacin rows in that file follow the OTHER (2nd-set)
    # layout and are byte-identical to the 2nd-set file, so they are skipped and
    # taken from df2 above.
    layout_1st = list(
        pd.read_csv(RAW / "mic" / "Summary_pathogens.csv", nrows=0).columns
    )[2:]
    df3 = pd.read_csv(
        RAW / "mic" / "MIC_Summary_Pathogens_with_lower_than_1uM_results.csv"
    )
    pids = df3[df3.columns[0]].astype(str).str.strip().str.strip('"')
    df3 = df3[df3.columns[1:]]
    assert len(df3.columns) == len(layout_1st), "sub-1 file: unexpected column count"
    df3.columns = layout_1st                      # relabel; shipped header is wrong

    ANTIBIOTIC_ROWS = {"Polymyxin B", "Levofloxacin"}
    seen_9310 = 0
    names_3rd, keep = [], []
    for pid in pids:
        if pid in ANTIBIOTIC_ROWS:
            names_3rd.append(None); keep.append(False); continue
        if pid == "sequence_9310":
            seen_9310 += 1
            names_3rd.append("Ω-DP-68" if seen_9310 == 1 else "Ω-MT-bZIP-4")
        else:
            names_3rd.append(resolve(pid, "", by_seqid, by_seq))
        keep.append(True)
    df3 = df3.assign(short_name=names_3rd).loc[keep]
    overlay = (
        df3.reindex(columns=["short_name"] + MIC_PATHOGENS)
        .drop_duplicates(subset=["short_name"], keep="first")
        .set_index("short_name")
    )

    base = combined.set_index("short_name")
    order = list(base.index) + [n for n in overlay.index if n not in set(base.index)]
    base = base.reindex(order)
    base.loc[overlay.index, MIC_PATHOGENS] = overlay[MIC_PATHOGENS]
    combined = (
        base.reindex(columns=MIC_PATHOGENS).rename_axis("short_name").reset_index()
    )

    # ── Censor NaN → ">64" ───────────────────────────────────────────────────
    for col in MIC_PATHOGENS:
        if col in combined.columns:
            combined[col] = combined[col].apply(censor_mic)

    combined.to_csv(OUT / "mic.csv", index=False)
    print(f"mic.csv:        {len(combined):>4} rows")


def build_disc(by_seqid, by_seq):
    """DiSC3-5 at-MIC: Peptide, Sequence, MaxRel, AUC."""
    df = pd.read_csv(
        RAW / "membrane_mechanism" / "disc" / "DiSC3-5_MaxRel_vs_AUC_atMIC.csv"
    )
    df["short_name"] = [
        resolve(p, str(s) if pd.notna(s) else "", by_seqid, by_seq)
        for p, s in zip(df["Peptide"], df["Sequence"])
    ]
    df[["short_name", "MaxRel", "AUC"]].to_csv(OUT / "disc.csv", index=False)
    print(f"disc.csv:       {len(df):>4} rows")


def _load_fc_concentrations(by_seqid, by_seq):
    """Map short_name -> tested fixed concentration (μM) for the MoA FC assays.

    Source: membrane_mechanism/tested_concentrations_FC.csv (per-peptide table).
    The concentration is family-anchored, so all analogs of a family share one
    value (e.g. cecropin/LG21/pa4/sarcotoxin = 0.5, bZIP = 2, Mammutin-1 = 16).
    Returns {} if the table is absent (concentration_uM then left blank).
    """
    path = RAW / "membrane_mechanism" / "tested_concentrations_FC.csv"
    if not path.exists():
        print("  WARNING: tested_concentrations_FC.csv missing — concentration_uM left blank")
        return {}
    t = pd.read_csv(path)
    conc = {}
    for _, r in t.iterrows():
        seq = str(r["Sequence"]) if pd.notna(r.get("Sequence")) else ""
        conc[resolve(str(r["Peptide"]), seq, by_seqid, by_seq)] = r["concentration_uM"]
    return conc


def build_disc_fc(by_seqid, by_seq):
    """DiSC3-5 fixed-concentration: blank id, MaxRel, AUC, per-peptide conc."""
    df = pd.read_csv(
        RAW / "membrane_mechanism" / "disc" / "DiSC3-5_MaxRel_vs_AUC_FC.csv"
    )
    df = df.rename(columns={df.columns[0]: "pid"})
    df["short_name"] = [resolve(p, "", by_seqid, by_seq) for p in df["pid"]]
    df = df[["short_name", "MaxRel", "AUC"]]
    df["concentration_uM"] = df["short_name"].map(_load_fc_concentrations(by_seqid, by_seq))
    df.to_csv(OUT / "disc_fc.csv", index=False)
    print(f"disc_fc.csv:    {len(df):>4} rows")


def build_npn(by_seqid, by_seq):
    """NPN at-MIC: Peptide, Sequence, MaxRel, AUC."""
    df = pd.read_csv(
        RAW / "membrane_mechanism" / "npn" / "NPN_MaxRel_vs_AUC_atMIC.csv"
    )
    df["short_name"] = [
        resolve(p, str(s) if pd.notna(s) else "", by_seqid, by_seq)
        for p, s in zip(df["Peptide"], df["Sequence"])
    ]
    df[["short_name", "MaxRel", "AUC"]].to_csv(OUT / "npn.csv", index=False)
    print(f"npn.csv:        {len(df):>4} rows")


def build_npn_fc(by_seqid, by_seq):
    """NPN fixed-concentration: blank id, MaxRel, AUC, per-peptide conc."""
    df = pd.read_csv(
        RAW / "membrane_mechanism" / "npn" / "NPN_MaxRel_vs_AUC _FC.csv"
    )
    df = df.rename(columns={df.columns[0]: "pid"})
    df["short_name"] = [resolve(p, "", by_seqid, by_seq) for p in df["pid"]]
    df = df[["short_name", "MaxRel", "AUC"]]
    df["concentration_uM"] = df["short_name"].map(_load_fc_concentrations(by_seqid, by_seq))
    df.to_csv(OUT / "npn_fc.csv", index=False)
    print(f"npn_fc.csv:     {len(df):>4} rows")


def build_hc50(by_seqid, by_seq):
    """
    Source: hc50/HC50_non-linear_regression_*.csv
    Format: blank, 'Red blood cells' header; one HC50 per row.
    Censoring: value > 128 → '>128'.
    bZIP group: sequence_9310 → Ω-MT-bZIP-4 (via BZIP_OVERRIDE).
    """
    frames = []
    for path in sorted((RAW / "hc50").glob("HC50_non-linear_regression_*.csv")):
        is_bzip = "bZIP" in path.name
        df = pd.read_csv(path, header=0)
        df = df.rename(columns={df.columns[0]: "pid", df.columns[1]: "HC50"})
        df["short_name"] = [
            resolve(p, "", by_seqid, by_seq,
                    overrides=BZIP_OVERRIDE if is_bzip else None)
            for p in df["pid"]
        ]
        df["HC50"] = df["HC50"].apply(censor_128)
        frames.append(df[["short_name", "HC50"]])

    result = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["short_name"])
    result.to_csv(OUT / "hc50.csv", index=False)
    print(f"hc50.csv:       {len(result):>4} rows")


def build_cc50(by_seqid, by_seq):
    """
    Source: cc50/OmegAMP_mean_all_replicates_CC50.csv
    Format: peptide IDs in row 0 (wide header), CC50 values in row 1.
    Censoring: value > 128 → '>128'.
    sequence_9310 collision: same as MIC/HC50 — 1st occurrence → Ω-DP-68,
    2nd occurrence → Ω-MT-bZIP-4 (confirmed by bZIP-family neighbours).
    """
    df = pd.read_csv(RAW / "cc50" / "OmegAMP_mean_all_replicates_CC50.csv", header=None)
    ids  = df.iloc[0].tolist()
    vals = df.iloc[1].tolist()

    seen_9310 = 0
    rows = []
    for pid, val in zip(ids, vals):
        pid_str = str(pid)
        if pid_str == "sequence_9310":
            seen_9310 += 1
            short = "Ω-DP-68" if seen_9310 == 1 else "Ω-MT-bZIP-4"
        else:
            short = resolve(pid_str, "", by_seqid, by_seq)
        rows.append({"short_name": short, "CC50": censor_128(val)})
    pd.DataFrame(rows).to_csv(OUT / "cc50.csv", index=False)
    print(f"cc50.csv:       {len(rows):>4} rows")


def build_dna_binding(by_seqid, by_seq):
    """
    Sources: dna_binding/delta*.txt  (tab-separated, wide format)
    All peptides in these files are bZIP analogs → apply BZIP_OVERRIDE.
    Output: long format (short_name, metric, concentration, value).

    deltaA / deltaLambdaNeg / deltaLambdaPos:
        rows = concentrations (12.5, 25.0, 50.0 μM)
        cols = peptides
    deltaAmax:
        rows = peptides, single value column (no concentration)
    """
    dna_dir = RAW / "dna_binding"

    def _resolve_bzip(pid):
        return resolve(str(pid), "", by_seqid, by_seq, overrides=BZIP_OVERRIDE)

    def parse_conc_file(path, metric):
        df = pd.read_table(path, index_col=0)
        df.index.name = "concentration"
        melted = (
            df.reset_index()
            .melt(id_vars="concentration", var_name="pid", value_name="value")
        )
        melted["short_name"] = melted["pid"].apply(_resolve_bzip)
        melted["metric"] = metric
        return melted[["short_name", "metric", "concentration", "value"]]

    def parse_amax_file(path):
        df = pd.read_table(path, header=0)
        df = df.rename(columns={df.columns[0]: "pid", df.columns[1]: "value"})
        df["short_name"] = df["pid"].apply(_resolve_bzip)
        df["metric"] = "deltaAmax"
        df["concentration"] = float("nan")
        return df[["short_name", "metric", "concentration", "value"]]

    result = pd.concat(
        [
            parse_conc_file(dna_dir / "deltaA_vs_DNA_all_concentrations.txt", "deltaA"),
            parse_amax_file(dna_dir / "deltaAmax_vs_DNA.txt"),
            parse_conc_file(dna_dir / "deltaLambdaNeg_vs_DNA_all_concentrations.txt", "deltaLambdaNeg"),
            parse_conc_file(dna_dir / "deltaLambdaPos_vs_DNA_all_concentrations.txt", "deltaLambdaPos"),
        ],
        ignore_index=True,
    )
    result.to_csv(OUT / "dna_binding.csv", index=False)
    print(f"dna_binding.csv:{len(result):>4} rows")


def build_lps_binding(by_seqid, by_seq):
    """
    Sources: lps_binding/{Family}_at_45min[_Ca].txt  (tab-separated)
    Header row: blank, 32, 32, 32, 8, 8, 8, 2, 2, 2  (3 replicates per concentration)
    Output: long format (short_name, family, concentration, replicate, calcium, BC_displacement).
    """
    lps_dir = RAW / "lps_binding"
    CONCS = [32, 32, 32, 8, 8, 8, 2, 2, 2]
    REPS  = [ 1,  2,  3, 1, 2, 3, 1, 2, 3]
    rows = []

    for family_key, family_name in LPS_FAMILIES.items():
        for calcium, suffix in [("no", "_at_45min.txt"), ("yes", "_at_45min_Ca.txt")]:
            path = lps_dir / f"{family_key}{suffix}"
            if not path.exists():
                print(f"  WARNING: {path.name} not found — skipping")
                continue

            with open(path, newline="") as fh:
                reader = csv_module.reader(fh, delimiter="\t")
                next(reader)  # skip header (concentrations row)
                for line in reader:
                    if not line or not line[0]:
                        continue
                    pid_raw = line[0].strip()
                    # Normalise template names; then try seq_id resolution
                    short = LPS_NAME_NORM.get(pid_raw) or resolve(pid_raw, "", by_seqid, by_seq)
                    vals = line[1:]
                    for val, conc, rep in zip(vals, CONCS, REPS):
                        try:
                            bc = float(val)
                        except (ValueError, TypeError):
                            continue  # skip empty / non-numeric cells
                        rows.append({
                            "short_name":    short,
                            "family":        family_name,
                            "concentration": conc,
                            "replicate":     rep,
                            "calcium":       calcium,
                            "BC_displacement": bc,
                        })

    pd.DataFrame(rows).to_csv(OUT / "lps_binding.csv", index=False)
    print(f"lps_binding.csv:{len(rows):>4} rows")


def build_bestsel(by_seqid, by_seq):
    """
    Primary source:  {Group}_BeStSel_{H2O|TFE|SDS|MeOH}_heatmap.csv
                     Four per-solvent files merged per group → 20-col rows.
    Secondary:       *_all_sec_str_all_solv_heatmap.csv  (pre-merged, 20 cols)

    Raw column order per fraction group: H2O, TFE/H2O, SDS/H2O, MeOH/H2O.
    NOTE: the previous bestsel.csv in the repo had TFE_H2O and MeOH_H2O
    labels swapped; this script assigns them correctly from the raw file order.
    """
    ss_dir = RAW / "secondary_structure"
    frames = []
    csv_groups: set[str] = set()

    # ── Per-solvent CSV groups (merge 4 solvents → 20 cols) ──────────────────
    for h2o_path in sorted(ss_dir.glob("*_BeStSel_H2O_heatmap.csv")):
        prefix = h2o_path.name.removesuffix("_BeStSel_H2O_heatmap.csv")
        paths = [ss_dir / f"{prefix}{suf}" for suf in _SS_SOLVENT_SUFFIXES]
        if not all(p.exists() for p in paths):
            missing = [p.name for p in paths if not p.exists()]
            print(f"  WARNING: {prefix}: missing solvent files {missing} — skipping group")
            continue

        overrides = BZIP_OVERRIDE if prefix == "bZIP" else None
        merged: pd.DataFrame | None = None
        for path, solvent in zip(paths, BESTSEL_SOLVENTS):
            df = pd.read_csv(path, header=0)
            df = df.rename(columns={df.columns[0]: "pid"})
            frac_cols = list(df.columns[1:6])
            rename_map = {c: f"{f}_{solvent}" for c, f in zip(frac_cols, BESTSEL_FRACTIONS)}
            df = df.rename(columns=rename_map)[["pid"] + [f"{f}_{solvent}" for f in BESTSEL_FRACTIONS]]
            merged = df if merged is None else merged.merge(df, on="pid", how="outer")

        assert merged is not None
        merged["short_name"] = [
            resolve(p, "", by_seqid, by_seq, overrides=overrides) for p in merged["pid"]
        ]
        frames.append(merged[["short_name"] + BESTSEL_COLS])
        csv_groups.add(prefix)
        print(f"  CSV group {prefix}: {len(merged)} rows")

    # ── All-in-one CSVs (20 cols already) ────────────────────────────────────
    for path in sorted(ss_dir.glob("*.csv")):
        if _SS_PER_SOLVENT_RE.search(path.name):
            continue  # handled above
        df = pd.read_csv(path, header=0)
        df = df.rename(columns={df.columns[0]: "pid"})
        data_cols = list(df.columns[1:])
        if len(data_cols) != 20:
            print(f"  WARNING: {path.name} has {len(data_cols)} data columns (expected 20) — skipping")
            continue
        df = df.rename(columns=dict(zip(data_cols, BESTSEL_COLS)))
        df["short_name"] = [resolve(p, "", by_seqid, by_seq) for p in df["pid"]]
        frames.append(df[["short_name"] + BESTSEL_COLS])

    if not frames:
        print("  WARNING: no BeStSel data found — skipping bestsel.csv")
        return

    result = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["short_name"])

    # Guard against silently building a truncated dataset (see EXPECTED_DENOVO_BESTSEL).
    ref = pd.read_csv(REF_TABLE, usecols=["short_name", "category"])
    ref["short_name"] = ref["short_name"].str.strip()
    names = result["short_name"].str.strip()
    n_denovo = int(ref.loc[ref["short_name"].isin(names), "category"].eq("de_novo").sum())
    if n_denovo < EXPECTED_DENOVO_BESTSEL:
        raise SystemExit(
            f"  ERROR: bestsel.csv has {n_denovo} de novo peptides "
            f"(expected {EXPECTED_DENOVO_BESTSEL}). Per-solvent BeStSel CSVs are "
            f"likely missing from {RAW / 'secondary_structure'} — refusing to "
            f"write a truncated dataset."
        )

    result.to_csv(OUT / "bestsel.csv", index=False)
    print(f"bestsel.csv:    {len(result):>4} rows  ({n_denovo} de novo)")


def build_proteolysis():
    """
    Source: proteolytic_resistance/Percentage_of_remaining_peptides.csv
    Format: 'Time (min)' column, then pairs of identically-named columns (2 replicates).
    Already uses short_names — no ID mapping needed.
    Output: long format (short_name, time_min, rep1, rep2, mean).
    """
    path = RAW / "proteolytic_resistance" / "Percentage_of_remaining_peptides.csv"

    with open(path, newline="") as fh:
        reader = csv_module.reader(fh)
        header = next(reader)
        data_rows = [row for row in reader if row]

    times = [float(row[0]) for row in data_rows]
    cols = header[1:]  # strip the "Time (min)" column

    rows = []
    i = 0
    while i < len(cols):
        name = cols[i]
        # Each peptide occupies exactly 2 consecutive columns (2 replicates)
        if i + 1 < len(cols) and cols[i + 1] == name:
            for t, row in zip(times, data_rows):
                r1 = float(row[i + 1])   # +1 offset: col 0 = Time
                r2 = float(row[i + 2])
                rows.append({
                    "short_name": name,
                    "time_min":   t,
                    "rep1":       round(r1, 1),
                    "rep2":       round(r2, 1),
                    "mean":       round((r1 + r2) / 2, 1),
                })
            i += 2
        else:
            i += 1  # orphan column — skip

    pd.DataFrame(rows).to_csv(OUT / "proteolysis.csv", index=False)
    print(f"proteolysis.csv:{len(rows):>4} rows")


# raw murine group label -> canonical short_name. bZIP is a motif design (MT),
# not an analog; the raw files label it AMT, so it is corrected here.
MURINE_NAMES = {
    "Control": "Control", "Polymyxin B": "Polymyxin B", "Levofloxacin": "Levofloxacin",
    "AT-BoCo1-5": "Ω-AT-BoCo1-5", "AT-BoCo1-9": "Ω-AT-BoCo1-9",
    "DP-52": "Ω-DP-52", "DP-19": "Ω-DP-19",
    "AMT-cecropin-1": "Ω-AMT-cecropin-1", "AMT-cecropin-4": "Ω-AMT-cecropin-4",
    "AMT-bZIP-8": "Ω-MT-bZIP-8", "AMT-pa4-1": "Ω-AMT-pa4-1",
}
MURINE = {
    "murine_skin.csv": {
        "cfu": [("CFU - Skin Scarification - A. baumannii ATCC19606 - Day 2.csv", 2),
                ("CFU - Skin Scarification - A. baumannii ATCC19606 - Day 4.csv", 4)],
        "weight": "Normalize of Mouse weight - Skin Scarification - A. baumannii ATCC19606.csv"},
    "murine_thigh.csv": {
        "cfu": [("CFU - Thigh Infection - A. baumannii ATCC19606 - Day 6.csv", 6),
                ("CFU - Thigh Infection - A. baumannii ATCC19606 - Day 8.csv", 8)],
        "weight": "Normalize of Mouse weight - Thigh Infection - A. baumannii ATCC19606.csv"},
}


def build_murine():
    """In vivo CFU (raw CFU g-1, one file per day) and body-weight (%) time-courses,
    one tidy file per model: short_name, assay (cfu|weight), day, replicate, value.
    """
    d = RAW / "murine"
    for out, spec in MURINE.items():
        rows = []
        for fname, day in spec["cfu"]:                       # CFU: one data row per file
            rr = list(csv_module.reader(open(d / fname)))
            hdr, vals = rr[0], rr[1]
            rep = {}
            for i in range(1, len(hdr)):
                v = vals[i].strip()
                if not v:
                    continue
                sn = MURINE_NAMES.get(hdr[i].strip(), hdr[i].strip())
                rep[sn] = rep.get(sn, 0) + 1
                rows.append({"short_name": sn, "assay": "cfu", "day": day,
                             "replicate": rep[sn], "value": float(v)})
        rr = list(csv_module.reader(open(d / spec["weight"])))  # weight: time-course
        hdr = rr[0]
        for r in rr[1:]:
            if not r or not r[0].strip():
                continue
            day = int(float(r[0]))
            rep = {}
            for i in range(1, len(hdr)):
                v = r[i].strip()
                if not v:
                    continue
                sn = MURINE_NAMES.get(hdr[i].strip(), hdr[i].strip())
                rep[sn] = rep.get(sn, 0) + 1
                rows.append({"short_name": sn, "assay": "weight", "day": day,
                             "replicate": rep[sn], "value": float(v)})
        pd.DataFrame(rows)[["short_name", "assay", "day", "replicate", "value"]].to_csv(OUT / out, index=False)
        print(f"{out}:{len(rows):>4} rows")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if not RAW.exists():
        print(f"ERROR: raw data directory not found at {RAW}", file=sys.stderr)
        print("Place Marcelo's raw files under data/raw/ and re-run.", file=sys.stderr)
        sys.exit(1)

    by_seqid, by_seq = load_ref()

    build_mic(by_seqid, by_seq)
    build_disc(by_seqid, by_seq)
    build_disc_fc(by_seqid, by_seq)
    build_npn(by_seqid, by_seq)
    build_npn_fc(by_seqid, by_seq)
    build_hc50(by_seqid, by_seq)
    build_cc50(by_seqid, by_seq)
    build_dna_binding(by_seqid, by_seq)
    build_lps_binding(by_seqid, by_seq)
    build_bestsel(by_seqid, by_seq)
    build_proteolysis()
    build_murine()
