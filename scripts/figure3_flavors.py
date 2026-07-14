"""
Loader for the figure-3 flavors data layout (data/figure3/flavors).

Mirrors the philosophy of Panels A/B (figure3_sweep): the OmegAMP flavors are
generated under the omegamp-inference API and read from ``data_root``, while the
external baselines — the fixed Prototypes set and the HydrAMP methods — are
taken from ``ref_root`` (the same path, post-consolidation). The split is kept
because the two use different file-naming conventions.

The new layout differs from the reference in three ways the reference loader
cannot handle:

  * predictions are ``sequence, MIC, MIC_unit`` (no ``Sequence_id`` column) and
    are deduplicated, so MIC is matched on the sequence *string*;
  * analog FASTA headers are positional (``>sequence_N``) with no prototype tag,
    so the analog->prototype mapping is recovered positionally, exactly like the
    sweep loader: prototypes are repeated R times in proto-major order, so
    ``prototype(i) = i // R`` with ``R = n_records // n_prototypes``;
  * de-novo similarity is named ``{stem}_vs_{pt}_prototypes.csv``.

Public entry point: :func:`load_flavors_new`, returning ``(mic_df, sim_df)`` with
the same schema as the notebook's original ``load_flavors_best`` —
``mic_df``: columns ``method, pt, mic``; ``sim_df``: columns ``method, pt, sim``.
"""
from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO


# ── new-data OmegAMP flavors (read from data_root) ───────────────────────────
_DENOVO = [  # (method_dir, stem, label)
    ("omegamp_du", "OmegAMP-U",     "OmegAMP-DU"),
    ("omegamp_dt", "OmegAMP-T",     "OmegAMP-DT"),
    ("omegamp_dp", "denovo_subset", "OmegAMP-DP"),
]
_ANALOG = [  # (method_dir, stem_tpl, label)
    ("omegamp_au",  "{pt}_analog_only",              "OmegAMP-AU"),
    ("omegamp_at",  "{pt}_analog_property",          "OmegAMP-AT"),
    ("omegamp_am",  "{pt}_analog_template",          "OmegAMP-AM"),
    ("omegamp_amt", "{pt}_analog_template_property", "OmegAMP-AMT"),
    ("omegamp_ap",  "{pt}_analog_subset",            "OmegAMP-AP"),
]

# ── reference baselines kept unchanged (read from ref_root) ──────────────────
_REF_NONANALOG_MIC = [  # (path_tpl, label, has_pt)
    ("prototypes/predictions/{pt}_prototypes-min.tsv", "Prototypes", True),
    ("hydramp_d/predictions/unconditional-min.tsv",    "HydrAMP-D",  False),
]
_REF_NONANALOG_SIM = [  # (pos_path, neg_path, label)
    ("prototypes/similarity/positive_prototypes_vs_curated-AMPs_samples.csv",
     "prototypes/similarity/negative_prototypes_vs_curated-Non-AMPs_samples.csv", "Prototypes"),
    ("hydramp_d/similarity/unconditional_vs_positive_prototypes.csv",
     "hydramp_d/similarity/unconditional_vs_negative_prototypes.csv", "HydrAMP-D"),
]
_REF_ANALOG = [  # (method_dir, stem_tpl, label, proto_regex)
    ("hydramp_a_t1",   "{pt}_1",   "HydrAMP-A τ=1",   r"prototype=([A-Z]+)"),
    ("hydramp_a_t2p5", "{pt}_2.5", "HydrAMP-A τ=2.5", r"prototype=([A-Z]+)"),
    ("hydramp_a_t5",   "{pt}_5",   "HydrAMP-A τ=5",   r"prototype=([A-Z]+)"),
]


def _n_unique(fasta):
    ids = set()
    for r in SeqIO.parse(str(fasta), "fasta"):
        ids.add(r.id)
    return len(ids)


def _load_denovo(F, mic_rows, sim_rows, pred_subdir="predictions"):
    """De-novo flavors: keep all predictions; MIC duplicated into both panels,
    similarity split by the vs-positive / vs-negative prototype file."""
    for mdir, stem, label in _DENOVO:
        pred = F / mdir / pred_subdir / f"{stem}.tsv"
        if pred.exists():
            mic = pd.read_csv(pred, sep="\t")["MIC"]
            for pt in ("positive", "negative"):
                for v in mic:
                    mic_rows.append({"method": label, "pt": pt, "mic": v})
        for pt in ("positive", "negative"):
            sim = F / mdir / "similarity" / f"{stem}_vs_{pt}_prototypes.csv"
            if not sim.exists():
                continue
            sdf = pd.read_csv(sim)
            for v in sdf.groupby("queryId")["normalized"].max().values:
                sim_rows.append({"method": label, "pt": pt, "sim": v})


def _load_analog(F, mic_rows, sim_rows, n_proto, pred_subdir="predictions"):
    """Analog flavors: lowest-MIC analog per prototype (positional mapping)."""
    for mdir, stem_tpl, label in _ANALOG:
        for pt in ("positive", "negative"):
            stem = stem_tpl.format(pt=pt)
            fasta = F / mdir / f"{stem}.fasta"
            pred = F / mdir / pred_subdir / f"{stem}.tsv"
            sim = F / mdir / "similarity" / f"{stem}_vs_{pt}_prototypes.csv"
            if not fasta.exists() or not pred.exists():
                continue

            records = list(SeqIO.parse(str(fasta), "fasta"))
            reps = max(1, len(records) // n_proto) if n_proto else 1
            mic_by_seq = pd.read_csv(pred, sep="\t").groupby("sequence")["MIC"].min().to_dict()
            sim_by_id = {}
            if sim.exists():
                sdf = pd.read_csv(sim)
                if "normalized" in sdf.columns and not sdf.empty:
                    sim_by_id = sdf.groupby("queryId")["normalized"].max().to_dict()

            # one champion (lowest MIC) per prototype group
            best = {}  # proto_idx -> (mic, seq_id)
            for i, r in enumerate(records):
                mic = mic_by_seq.get(str(r.seq))
                if mic is None:
                    continue
                pidx = i // reps
                if pidx not in best or mic < best[pidx][0]:
                    best[pidx] = (mic, r.id)
            for mic, sid in best.values():
                mic_rows.append({"method": label, "pt": pt, "mic": mic})
                if sid in sim_by_id:
                    sim_rows.append({"method": label, "pt": pt, "sim": sim_by_id[sid]})


def _load_ref_baselines(R, mic_rows, sim_rows):
    """Prototypes + HydrAMP baselines (APEX).

    MIC from the fixed reference ``*-min.tsv`` files; similarity from the same
    reference CSVs. HydrAMP-A keeps the per-prototype champion (lowest MIC).
    """
    # ── non-analog MIC (Prototypes, HydrAMP-D) ────────────────────────────────
    for path_tpl, label, has_pt in _REF_NONANALOG_MIC:
        for pt in (["positive", "negative"] if has_pt else ["positive"]):
            fp = R / (path_tpl.format(pt=pt) if has_pt else path_tpl)
            if not fp.exists():
                continue
            pdf = pd.read_csv(fp, sep="\t", index_col=0)
            for t in ([pt] if has_pt else ["positive", "negative"]):
                for v in pdf["MIC"]:
                    mic_rows.append({"method": label, "pt": t, "mic": v})

    # ── non-analog similarity ─────────────────────────────────────────────────
    for pos_p, neg_p, label in _REF_NONANALOG_SIM:
        for pt, path in [("positive", pos_p), ("negative", neg_p)]:
            fp = R / path
            if not fp.exists():
                continue
            sdf = pd.read_csv(fp)
            for v in sdf.groupby("queryId")["normalized"].max().values:
                sim_rows.append({"method": label, "pt": pt, "sim": v})

    # ── HydrAMP-A analog MIC (champion per prototype) + champion similarity ────
    for mdir, stem_tpl, label, proto_regex in _REF_ANALOG:
        best_ids = {}
        for pt in ("positive", "negative"):
            stem = stem_tpl.format(pt=pt)
            tsv = R / mdir / "predictions" / f"{stem}-min.tsv"
            if not tsv.exists():
                continue
            pdf = pd.read_csv(tsv, sep="\t", index_col=0)
            pdf["proto_key"] = pdf["Sequence_id"].apply(
                lambda x: m.group(1) if (m := re.search(proto_regex, str(x))) else None)
            pdf = pdf.dropna(subset=["proto_key"])
            bst = pdf.loc[pdf.groupby("proto_key")["MIC"].idxmin()]
            for _, row in bst.iterrows():
                mic_rows.append({"method": label, "pt": pt, "mic": row["MIC"]})
            best_ids[pt] = set(bst["Sequence_id"].values)
        for pt in ("positive", "negative"):
            stem = stem_tpl.format(pt=pt)
            sim_csv = R / mdir / "similarity" / f"{stem}_vs_{pt}_prototypes.csv"
            if not sim_csv.exists():
                continue
            sdf = pd.read_csv(sim_csv)
            sdf = sdf[sdf["queryId"].isin(best_ids.get(pt, set()))]
            for v in sdf.groupby("queryId")["normalized"].max().values:
                sim_rows.append({"method": label, "pt": pt, "sim": v})


def load_flavors_new(data_root, ref_root, pos_prototypes, neg_prototypes):
    """Return ``(mic_df, sim_df)`` for the new flavors run (APEX-scored).

    OmegAMP flavors come from ``data_root``; Prototypes + HydrAMP baselines come
    from ``ref_root`` unchanged. ``pos_prototypes`` / ``neg_prototypes`` are the
    prototype FASTAs, used only to count unique prototypes for positional
    analog->prototype mapping (positive and negative sets are the same size).
    """
    F = Path(data_root) / "figure3/flavors"
    R = Path(ref_root) / "figure3/flavors"
    n_proto = _n_unique(pos_prototypes)

    mic_rows, sim_rows = [], []
    _load_denovo(F, mic_rows, sim_rows)
    _load_analog(F, mic_rows, sim_rows, n_proto)
    _load_ref_baselines(R, mic_rows, sim_rows)
    return pd.DataFrame(mic_rows), pd.DataFrame(sim_rows)
