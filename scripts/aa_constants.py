"""
Amino acid physicochemical constants shared across all figure scripts and notebooks.
"""

WATER = 18.015  # molecular weight of water (Da), subtracted per peptide bond

AA_MW = {
    'A':  89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'E': 147.13, 'Q': 146.15, 'G':  75.03, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
}

# Eisenberg consensus hydrophobicity scale
EISENBERG = {
    'A':  0.620, 'R': -2.530, 'N': -0.780, 'D': -0.900, 'C':  0.290,
    'Q': -0.850, 'E': -0.740, 'G':  0.480, 'H': -0.400, 'I':  1.380,
    'L':  1.060, 'K': -1.500, 'M':  0.640, 'F':  1.190, 'P':  0.120,
    'S': -0.180, 'T': -0.050, 'W':  0.810, 'Y':  0.260, 'V':  1.080,
}

# AGGRESCAN aggregation propensity scale (de Groot et al., 2006)
AGGRESCAN = {
    'I':  1.822, 'F':  1.754, 'V':  1.594, 'L':  1.380, 'Y':  1.159,
    'W':  1.037, 'M':  0.910, 'C':  0.604, 'A': -0.036, 'T': -0.159,
    'S': -0.294, 'P': -0.334, 'G': -0.535, 'K': -0.931, 'H': -1.033,
    'Q': -1.231, 'R': -1.240, 'N': -1.302, 'E': -1.412, 'D': -1.836,
}

# Net charge contribution at pH 7 (simplified: K/R +1, D/E -1)
AA_CHARGE = {'K': 1, 'R': 1, 'D': -1, 'E': -1}


def peptide_mw(seq: str) -> float:
    """Monoisotopic MW of a peptide (Da): sum of residue MWs minus water per bond."""
    s = str(seq).upper()
    return sum(AA_MW.get(aa, 0) for aa in s) - (len(s) - 1) * WATER
