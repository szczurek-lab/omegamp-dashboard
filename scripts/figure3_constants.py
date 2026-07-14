"""
Shared constants for Figure 3 (analog generation sweep analysis).

Import in notebooks with:
    sys.path.insert(0, os.path.abspath('..'))
    from scripts.figure3_constants import ...
"""
import numpy as np

# ── Sweep parameter grid ──────────────────────────────────────────────────────
TAU_VALUES   = [0, 50, 100, 150, 200, 300, 500, 700]
SIGMA_VALUES = [0.0, 0.25, 0.5, 0.75, 1.0]

# Subset plotted in legends and line plots (sparse enough to read)
TAU_PLOT = [0, 100, 200, 300, 500, 700]

# ── Colour maps ───────────────────────────────────────────────────────────────
# Viridis sequence: cold = conservative (low τ), warm = diverse (high τ)
TAU_COLORS = {
      0: '#440154',
    100: '#365c8d',
    200: '#1fa187',
    300: '#4ac16d',
    500: '#9fda3a',
    700: '#fde725',
}

# Sequential red: light = low noise (σ=0), dark = high noise (σ=1)
SIGMA_COLORS = {
    0.00: '#fee5d9',
    0.25: '#fcae91',
    0.50: '#fb6a4a',
    0.75: '#de2d26',
    1.00: '#a50f15',
}

# Scatter marker area proportional to τ magnitude
TAU_SIZES = {0: 6, 50: 10, 100: 14, 150: 18, 200: 22, 300: 28, 500: 36, 700: 44}

# ── Activity / similarity thresholds ─────────────────────────────────────────
# MIC_THRESH: log₂(32 µM) — upper boundary of the "active" green zone in Panel B.
# SIM_THRESH: minimum prototype-similarity required to stay in the green zone.
LOG32      = float(np.log2(32))  # ≈ 5.0
MIC_THRESH = LOG32
SIM_THRESH = 0.6

# ── Count-bar normalisation ───────────────────────────────────────────────────
# De novo methods produce ~50 k sequences; analog methods produce ~500.
# Each group is normalised independently so bars are comparable within group.
DENOVO_MAX = 50_000
ANALOG_MAX = 500

# ── Visual axis-break geometry (data coordinates) ────────────────────────────
# Covers the gap between analog methods (x ≤ 10) and de novo methods (x ≥ 13)
# in Panels C and D.
BREAK_X0 = 10.6
BREAK_X1 = 12.4

# ── Method display order, grouping, and colours ───────────────────────────────
# '_sep_' entries are rendered as dotted separator lines (not a real method).
METHOD_ORDER = [
    'Prototypes', '_sep_',
    'OmegAMP-AU', 'OmegAMP-AT', 'OmegAMP-AP', 'OmegAMP-AM', 'OmegAMP-AMT', '_sep_',
    'HydrAMP-A τ=1', 'HydrAMP-A τ=2.5', 'HydrAMP-A τ=5', '_sep_', '_sep_',
    'OmegAMP-DU', 'OmegAMP-DT', 'OmegAMP-DP', '_sep_',
    'HydrAMP-D',
]
SEP_POSITIONS = [1, 7, 11, 12, 16]

# Per-method violin/bar fill colours
GRP_COLORS = {
    'OmegAMP-DU':      '#F4A89A',
    'OmegAMP-DT':      '#E05C4B',
    'OmegAMP-DP':      '#9E2A2B',
    'HydrAMP-D':       '#5BA4A4',
    'Prototypes':      '#808080',
    'OmegAMP-AU':      '#D4B0B0',
    'OmegAMP-AT':      '#C9A0A0',
    'OmegAMP-AM':      '#A35D5D',
    'OmegAMP-AMT':     '#7B2D3D',
    'OmegAMP-AP':      '#5C1A2A',
    'HydrAMP-A τ=1':   '#7FCAC3',
    'HydrAMP-A τ=2.5': '#3A9E95',
    'HydrAMP-A τ=5':   '#1A5F5A',
}

# Methods that generate sequences without a prototype; used for count-bar
# normalisation (DENOVO_MAX vs ANALOG_MAX).
DENOVO_METHODS = frozenset({'OmegAMP-DU', 'OmegAMP-DT', 'OmegAMP-DP', 'HydrAMP-D'})
