# amp_similarity

AMP-specific sequence similarity scoring using a custom substitution matrix derived from antimicrobial peptide alignments.

## Usage

```python
from amp_similarity import normalized_alignment_score

score = normalized_alignment_score("KWKLFKKIEKVGRNIR", "KWRLFRKIEKVGRNIR")
# Returns 0-100 (percentage of self-alignment score)
```

The score is computed as:

    score(query, target) / score(query, query) * 100

using Biopython's `PairwiseAligner` in local alignment mode with the AMP substitution matrix (`amp_matrix.txt`), gap open penalty -10, gap extend penalty -1.

## Dependencies

- biopython
- numpy

## The AMP substitution matrix

`amp_matrix.txt` is a 20x20 log-odds matrix derived from curated AMP sequence alignments using a BLOSUM-like procedure:

1. AMP sequences are clustered with CD-HIT
2. Each cluster is aligned with MAFFT using a preliminary AMP matrix
3. Ungapped blocks are extracted from each alignment
4. Amino acid substitution pairs are counted across all blocks (BLOSUM counting at 100% threshold)
5. Counts are converted to log-odds scores

The full derivation pipeline is in `derive_matrix/`. To reproduce:

```bash
# Requires: cd-hit, mafft, pysum, biopython
cd derive_matrix
python run_pipeline.py -f <amp_fasta> -c <cdhit_clstr> -o output
# Then run combine_counters_and_derive_matrix.ipynb to produce the final matrix
```

Source: AMP_Metrics_Evaluation (CITE).
