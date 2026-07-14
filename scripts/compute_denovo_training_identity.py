#!/usr/bin/env python3
"""
Compute max length-matched sequence identity between each de novo peptide
and the training set (IsAMP == 1).

Identity is defined as:
    identical_positions / max(len_query, len_target)
using a global alignment (match=1, mismatch=0, gap=-2).
Only training sequences within ±40 % of the query length are considered.

Output
------
data/denovo_training_identity.csv
    short_name, sequence, length, best_id_lenmatched, best_hit_sequence
"""

import os
import sys
import time
import pandas as pd
from Bio.Align import PairwiseAligner

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)

DATA = os.path.join(_ROOT, 'data')
TRAIN_CSV = os.path.join(DATA, 'figure2', 'training', 'generative-model-dataset.csv')
REF_CSV   = os.path.join(DATA, 'omegamp_reference_table.csv')
OUT_CSV   = os.path.join(DATA, 'denovo_training_identity.csv')

_aligner = PairwiseAligner()
_aligner.mode = 'global'
_aligner.match_score = 1
_aligner.mismatch_score = 0
_aligner.open_gap_score = -2
_aligner.extend_gap_score = -0.5


def identity(query, target):
    score = _aligner.score(query, target)
    return score / max(len(query), len(target))


def main():
    print('Loading reference table …')
    ref = pd.read_csv(REF_CSV)
    denovo = ref[ref['category'] == 'de_novo'][['short_name', 'sequence']].dropna()
    denovo['sequence'] = denovo['sequence'].str.strip()
    print(f'  {len(denovo)} de novo peptides')

    print(f'Loading training AMPs … ({TRAIN_CSV})')
    t0 = time.time()
    df = pd.read_csv(TRAIN_CSV)
    train_seqs = df[df['IsAMP'] == 1]['Sequence'].dropna().tolist()
    print(f'  {len(train_seqs):,} AMP sequences loaded in {time.time()-t0:.1f}s')

    rows = []
    for _, row in denovo.iterrows():
        name = row['short_name']
        seq  = row['sequence']
        qlen = len(seq)
        lo, hi = qlen * 0.6, qlen * 1.4

        best_id, best_seq = 0.0, ''
        for t_seq in train_seqs:
            if not (lo <= len(t_seq) <= hi):
                continue
            sid = identity(seq, t_seq)
            if sid > best_id:
                best_id, best_seq = sid, t_seq

        rows.append({
            'short_name':        name,
            'sequence':          seq,
            'length':            qlen,
            'best_id_lenmatched': round(best_id, 4),
            'best_hit_sequence': best_seq,
        })
        print(f'  [{name}]  {best_id*100:.1f}%  {best_seq}')

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f'\nSaved → {OUT_CSV}')
    print(out[['short_name', 'best_id_lenmatched']].to_string(index=False))


if __name__ == '__main__':
    main()