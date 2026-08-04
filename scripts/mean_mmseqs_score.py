"""
Compute mean normalized MMseqs2 bitscore between a query FASTA and a database FASTA.

Core function
-------------
compute_mean_normalized_bitscore(query_fasta, database_fasta, output_results_csv=None)
    Returns (mean_normalized_score, coverage).
    If output_results_csv is given, writes a per-query CSV with columns:
        queryId, targetId, bit_score, self_targetId, self_bit_score, normalized

Temp files (result_*.tsv, fasta_filtered.fasta, tmp/) are written to the current
working directory, so call from a dedicated or temporary directory.

CLI usage
---------
    python scripts/mean_mmseqs_score.py --query <fasta> --database <fasta>
"""

import subprocess
import tempfile
import pandas as pd
import os
import argparse
import shutil
from pathlib import Path
from Bio import SeqIO

import shutil as _shutil
import os as _os

def _find_mmseqs():
    # 1. explicit override via env var
    if _os.environ.get("MMSEQS_BIN"):
        return _os.environ["MMSEQS_BIN"]
    # 2. same directory as this Python interpreter (conda env bin/)
    py_bin = _os.path.dirname(_os.__file__)          # e.g. .../lib/python3.11
    conda_bin = _os.path.join(py_bin, "../../bin")   # .../bin
    candidate = _os.path.normpath(_os.path.join(conda_bin, "mmseqs"))
    if _os.path.isfile(candidate):
        return candidate
    # 3. PATH
    found = _shutil.which("mmseqs")
    if found:
        return found
    raise FileNotFoundError(
        "mmseqs binary not found. Install it (e.g. conda install -c bioconda mmseqs2) "
        "or set the MMSEQS_BIN environment variable."
    )

# Resolved lazily: notebooks that only read cached similarity CSVs (mmseqs_sim below)
# must import fine on a host without the binary; the error still fires on real use.
PROFILE = "easy-search"
MMSEQS_THREADS = int(_os.environ.get("MMSEQS_THREADS", "2"))


def get_num_sequences(fasta_path):
    return sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))


def safe_load_tsv(path):
    if os.path.getsize(path) == 0:
        return pd.DataFrame(columns=[0, 11])
    df = pd.read_csv(path, sep="\t", header=None)
    df[0] = df[0].astype(str)  # match FASTA IDs which are always strings
    return df


def run_mmseq(query_fasta, database, i, add_self_matches=False):
    cmd = [
        _find_mmseqs(),
        PROFILE,
        query_fasta,
        database,
        f"result_{i}.tsv",
        "tmp",
        "-v", "0",
        "-e", "1000",
        "--comp-bias-corr", "0",   # disable composition bias correction for short peptides
        "--threads", str(MMSEQS_THREADS),
    ]
    if add_self_matches:
        # Self-alignment: default prefilter is fine (every seq shares k-mers with itself)
        cmd += ["--add-self-matches", "--max-seqs", "1"]
    else:
        # DB search: bypass k-mer prefilter so short peptides (<20 aa) are not dropped
        cmd += ["--prefilter-mode", "2"]
    subprocess.run(cmd, check=True)


def compute_mean_normalized_bitscore(
    query_fasta: str,
    database_fasta: str,
    output_results_csv: str | None = None,
) -> tuple[float, float]:
    i = 0

    query_records = list(SeqIO.parse(query_fasta, "fasta"))
    db_records = list(SeqIO.parse(database_fasta, "fasta"))
    db_seqs = {str(r.seq) for r in db_records}

    total_queries = len(query_records)
    query_ids = [r.id for r in query_records]

    shared_records  = [r for r in query_records if str(r.seq) in db_seqs]
    filtered_records = [r for r in query_records if str(r.seq) not in db_seqs]
    shared_ids = {r.id for r in shared_records}

    results = pd.DataFrame({
        "queryId":       query_ids,
        "targetId":      None,
        "bit_score":     None,
        "self_targetId": None,
        "self_bit_score": None,
        "normalized":    0.0,
    }).set_index("queryId")

    results.loc[list(shared_ids), "normalized"]    = 1.0
    results.loc[list(shared_ids), "targetId"]      = list(shared_ids)
    results.loc[list(shared_ids), "self_targetId"] = list(shared_ids)

    if not filtered_records:
        if output_results_csv:
            results.reset_index().to_csv(output_results_csv, index=False)
        return 1.0, 1.0

    filtered_fasta = "fasta_filtered.fasta"
    SeqIO.write(filtered_records, filtered_fasta, "fasta")

    run_mmseq(filtered_fasta, database_fasta, i, add_self_matches=False)
    df_db = safe_load_tsv(f"result_{i}.tsv")

    if df_db.empty:
        mean_norm = results["normalized"].mean()
        coverage  = len(shared_ids) / total_queries
        if output_results_csv:
            results.reset_index().to_csv(output_results_csv, index=False)
        return mean_norm, coverage

    idx = df_db.groupby(0)[df_db.columns[-1]].idxmax()
    best_db_hits = df_db.loc[idx, [0, 1, df_db.columns[-1]]]
    best_db_hits.columns = ["queryId", "targetId", "bit_score"]
    i += 1

    run_mmseq(filtered_fasta, filtered_fasta, i, add_self_matches=True)
    df_self = safe_load_tsv(f"result_{i}.tsv")

    if df_self.empty:
        mean_norm = results["normalized"].mean()
        coverage  = len(shared_ids) / total_queries
        if output_results_csv:
            results.reset_index().to_csv(output_results_csv, index=False)
        return mean_norm, coverage

    idx = df_self.groupby(0)[df_self.columns[-1]].idxmax()
    best_self_hits = df_self.loc[idx, [0, 1, df_self.columns[-1]]]
    best_self_hits.columns = ["queryId", "self_targetId", "self_bit_score"]

    merged = best_db_hits.merge(best_self_hits, on="queryId", how="inner")
    merged["normalized"] = merged["bit_score"] / merged["self_bit_score"]

    if (merged["normalized"] > 1).any():
        print("Warning: some normalized bitscores > 1")

    results.update(merged.set_index("queryId"))

    mean_norm = results["normalized"].mean()
    coverage  = (results["normalized"] > 0).mean()
    assert len(results) == total_queries

    if output_results_csv:
        results.reset_index().to_csv(output_results_csv, index=False)

    return mean_norm, coverage


def mmseqs_sim(query_fasta, database_fasta, cache_csv) -> dict:
    """Return {queryId: max_normalized_bitscore}, computing and caching if needed.

    Runs in a temporary directory so MMseqs2 scratch files don't pollute CWD.
    Results are written to cache_csv and reused on subsequent calls.
    """
    cache_csv = Path(cache_csv).resolve()
    if not cache_csv.exists() or cache_csv.stat().st_size == 0:
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        orig = _os.getcwd()
        with tempfile.TemporaryDirectory(prefix='mmseqs_') as tmpdir:
            try:
                _os.chdir(tmpdir)
                compute_mean_normalized_bitscore(
                    str(Path(query_fasta).resolve()),
                    str(Path(database_fasta).resolve()),
                    str(cache_csv),
                )
            finally:
                _os.chdir(orig)
    df = pd.read_csv(cache_csv)
    return df.groupby('queryId')['normalized'].max().to_dict()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute mean normalized MMseqs2 bitscore."
    )
    parser.add_argument("--query",    required=True, help="Query FASTA")
    parser.add_argument("--database", required=True, help="Database FASTA")
    parser.add_argument("--output",   default=None,  help="Output CSV (optional)")
    args = parser.parse_args()

    mean, cov = compute_mean_normalized_bitscore(args.query, args.database, args.output)
    print(f"Mean normalized bitscore: {mean:.4f}  coverage: {cov:.4f}")
