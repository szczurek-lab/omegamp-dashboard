"""Load the AMP-derived substitution matrix."""

from pathlib import Path

import numpy as np
from Bio.Align import substitution_matrices

_MATRIX_PATH = Path(__file__).parent / "amp_matrix.txt"


def load_amp_matrix():
    """Load the 20x20 AMP substitution matrix from amp_matrix.txt.

    Returns a Bio.Align.substitution_matrices.Array suitable for
    use with Bio.Align.PairwiseAligner.
    """
    lines = _MATRIX_PATH.read_text().strip().split("\n")
    aas = lines[0].split()
    data = np.array([[float(v) for v in row.split()[1:]] for row in lines[1:]])
    return substitution_matrices.Array(alphabet="".join(aas), dims=2, data=data)
