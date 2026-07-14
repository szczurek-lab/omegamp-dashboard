"""
Loader for the figure-3 sweep data layout (data/figure3/sweep).

This omegamp-inference layout has header-less FASTAs and deassociated
predictions, so the analog->prototype mapping is recovered positionally:

  sequences/tau{τ}_sigma{σ}.fasta          headers are just ">sequence_N"
  predictions/tau{τ}_sigma{σ}.tsv          columns: sequence, MIC, MIC_unit
                                           (already APEX-min, deduplicated)
  similarity/tau{τ}_sigma{σ}_vs_{pt}_prototypes.csv
                                           queryId = sequence_N, normalized = sim
  sequences/prototypes-repeated.fasta      the conditioning prototypes, each
                                           repeated R times, in proto-major order

Because the generated FASTA is in proto-major order (R consecutive analogs per
prototype) and ``prototypes-repeated.fasta`` lists the prototypes in that same
order, we recover the analog→prototype mapping positionally:

    n_prototypes = number of *unique* prototype ids
    reps         = n_records // n_prototypes        (10 for full cells, 2 for
                                                     the deduplicated τ=0 cells)
    prototype(i) = unique_prototypes[i // reps]

This was verified against the deterministic τ=0 cells, where every analog equals
its mapped prototype (1000/1000 and 5000/5000).

The public entry point is :func:`load_sweep_cell`, returning one row per
generated sequence with its prototype, predicted MIC and similarity.
"""
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def _tau_stem(tau):
    """Reference uses integer τ (0..700); the new files use the fraction."""
    return "0" if tau == 0 else str(tau / 1000.0)


def _ordered_unique_prototypes(proto_fasta):
    """Prototype (id, seq) pairs in first-seen order, deduplicated by id."""
    seen, uniq = set(), []
    for r in SeqIO.parse(str(proto_fasta), "fasta"):
        if r.id not in seen:
            seen.add(r.id)
            uniq.append((r.id, str(r.seq)))
    return uniq


def load_sweep_cell(data_root, set_name, tau, sigma, pred_subdir="predictions"):
    """Return a per-sequence DataFrame for one (τ, σ) sweep cell, or ``None``
    if the cell's FASTA/predictions are missing.

    Columns: seq_id, seq, proto_idx, proto_id, proto_seq, mic, sim.
    ``mic`` / ``sim`` are NaN where a sequence has no prediction / similarity
    entry (predictions are deduplicated, so MIC is matched on the sequence
    string; similarity is matched on the FASTA id).

    MIC comes from the APEX ``predictions/`` directory.
    """
    base = Path(data_root) / "figure3/sweep" / set_name
    stem = f"tau{_tau_stem(tau)}_sigma{sigma}"
    fasta = base / "sequences" / f"{stem}.fasta"
    pred = base / pred_subdir / f"{stem}.tsv"
    sim = base / "similarity" / f"{stem}_vs_{set_name}_prototypes.csv"

    if not fasta.exists() or not pred.exists():
        return None

    records = list(SeqIO.parse(str(fasta), "fasta"))
    uniq = _ordered_unique_prototypes(base / "sequences" / "prototypes-repeated.fasta")
    n_proto = len(uniq)
    reps = max(1, len(records) // n_proto) if n_proto else 1

    pdf = pd.read_csv(pred, sep="\t")
    mic_by_seq = pdf.groupby("sequence")["MIC"].min().to_dict()

    sim_by_id = {}
    if sim.exists():
        sdf = pd.read_csv(sim)
        if "normalized" in sdf.columns and not sdf.empty:
            sim_by_id = sdf.groupby("queryId")["normalized"].max().to_dict()

    rows = []
    for i, r in enumerate(records):
        pidx = i // reps
        pid, pseq = uniq[pidx] if pidx < n_proto else (None, None)
        seq = str(r.seq)
        rows.append({
            "seq_id": r.id,
            "seq": seq,
            "proto_idx": pidx,
            "proto_id": pid,
            "proto_seq": pseq,
            "mic": mic_by_seq.get(seq),
            "sim": sim_by_id.get(r.id),
        })
    return pd.DataFrame(rows)


def best_per_prototype(cell_df):
    """Lowest-MIC analog per prototype (drops sequences without a MIC)."""
    valid = cell_df.dropna(subset=["mic"])
    if valid.empty:
        return valid
    return valid.loc[valid.groupby("proto_idx")["mic"].idxmin()]