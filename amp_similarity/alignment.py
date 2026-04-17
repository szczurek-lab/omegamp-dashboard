"""AMP-specific pairwise alignment using a custom substitution matrix.

The substitution matrix was derived from curated AMP sequence alignments
using a BLOSUM-like procedure (see derive_matrix/ for the full pipeline).
"""

from Bio.Align import PairwiseAligner

from .matrix import load_amp_matrix

_aligner = None


def _get_aligner():
    global _aligner
    if _aligner is None:
        _aligner = PairwiseAligner()
        _aligner.substitution_matrix = load_amp_matrix()
        _aligner.open_gap_score = -10
        _aligner.extend_gap_score = -1
        _aligner.mode = "local"
    return _aligner


def normalized_alignment_score(query, target):
    """Normalized local alignment score using the AMP substitution matrix.

    Score = aligner.score(query, target) / aligner.score(query, query) * 100

    Parameters
    ----------
    query : str
        Query sequence (typically the generated analog).
    target : str
        Target sequence (typically the prototype).

    Returns
    -------
    float
        Normalized score as percentage (0--100).
    """
    if not query or not target:
        return 0.0
    aligner = _get_aligner()
    self_score = aligner.score(query, query)
    if self_score == 0:
        return 0.0
    return aligner.score(query, target) / self_score * 100
