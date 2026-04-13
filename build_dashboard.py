#!/usr/bin/env python3
"""
Build OmegAMP interactive dashboard.

Reads processed CSV files and a template to produce a self-contained HTML
dashboard with embedded JSON data. Requires: Python 3, modlAMP.

Usage:
    pip install modlamp
    python build_dashboard.py --data-dir ./data --output dashboard.html

Required input files (under --data-dir):
    omegamp_reference_table.csv  - peptide metadata
    mic.csv                      - MIC values (20 strain columns)
    cc50.csv                     - cytotoxicity CC50
    hc50.csv                     - hemolysis HC50

Optional input files (loaded if present):
    disc.csv           - DiSC3(5) membrane depolarization (MaxRel, AUC)
    npn.csv            - NPN outer membrane permeabilization (MaxRel, AUC)
    bestsel.csv        - BeStSel secondary structure (5 fracs × 4 solvents)
    lps_binding.csv    - LPS neutralization (BC displacement, ± Ca²⁺)
    dna_binding.csv    - DNA binding (ΔA at multiple concentrations)

Help files (under --help-dir or <data-dir>/help/):
    *.md               - Markdown files rendered as collapsible info panels.
                         File stem becomes the key (e.g. lps-binding.md -> "lps-binding").
                         Matched to help buttons via data-help attribute in the template.

Template:
    dashboard_template.html — must contain a single __DATA__ placeholder.
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import Counter, defaultdict

from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor


# ══════════════════════════════════════════════════════════════════════════
# Physicochemical descriptors (via modlAMP)
# ══════════════════════════════════════════════════════════════════════════

def compute_descriptors(seq):
    """
    Compute physicochemical descriptors for a peptide sequence using modlAMP.
    Returns dict with short keys for compact JSON, or all-None if seq is empty.

    Keys:
        q  — net charge at pH 7.4 (Henderson-Hasselbalch)
        h  — mean hydrophobicity (Eisenberg consensus)
        hm — hydrophobic moment (α-helix, 100° window)
        ar — aromaticity (fraction F/W/Y)
        ai — aliphatic index (Ikai 1980)
        fc — fraction charged residues (D/E/K/R/H)
        pI — isoelectric point
    """
    if not seq:
        return {k: None for k in ('q', 'h', 'hm', 'ar', 'ai', 'fc', 'pI')}

    g = GlobalDescriptor([seq])

    g.calculate_charge(ph=7.4, amide=False)
    charge = round(g.descriptor[0][0], 2)

    g.isoelectric_point(amide=False)
    pI = round(g.descriptor[0][0], 2)

    g.aliphatic_index()
    ai = round(g.descriptor[0][0], 2)

    g.aromaticity()
    ar = round(g.descriptor[0][0], 4)

    # Fraction charged (D, E, K, R, H)
    fc = round(sum(1 for a in seq if a in 'DEKRH') / len(seq), 4)

    # Hydrophobicity and moment via Eisenberg scale
    p = PeptideDescriptor([seq], 'eisenberg')
    p.calculate_global()
    h = round(p.descriptor[0][0], 4)

    p2 = PeptideDescriptor([seq], 'eisenberg')
    p2.calculate_moment(window=len(seq), angle=100)
    hm = round(p2.descriptor[0][0], 4)

    return {'q': charge, 'h': h, 'hm': hm, 'ar': ar, 'ai': ai, 'fc': fc, 'pI': pI}


def geometric_mean(values):
    """Geometric mean of non-None positive values. Returns None if empty."""
    v = [x for x in values if x is not None and x > 0]
    return round(math.exp(sum(math.log(x) for x in v) / len(v)), 3) if v else None


# ══════════════════════════════════════════════════════════════════════════
# Strain metadata — 20 strains matching mic.csv column order
# ══════════════════════════════════════════════════════════════════════════
# s = display name, g = Gram stain (+ or −), m = multidrug-resistant
#
# Gram assignments follow standard microbiology.
# MDR designations match ATCC catalog entries:
#   CRAB  = carbapenem-resistant A. baumannii (BAA-1605)
#   CRE   = carbapenem-resistant Enterobacteriaceae (AIC222, BAA-3170)
#   ESBL  = extended-spectrum β-lactamase K. pneumoniae (BAA-2342)
#   FQR   = fluoroquinolone-resistant P. aeruginosa (BAA-3197)
#   MRSA  = methicillin-resistant S. aureus (BAA-1556)
#   VRE   = vancomycin-resistant Enterococcus (700802, 700221)

STRAIN_INFO = [
    {"s": "A. baumannii ATCC 19606",       "g": "-", "m": False},
    {"s": "A. baumannii BAA-1605 (CRAB)",   "g": "-", "m": True},
    {"s": "E. cloacae ATCC 13047",          "g": "-", "m": False},
    {"s": "E. coli ATCC 11775",             "g": "-", "m": False},
    {"s": "E. coli AIC221",                 "g": "-", "m": False},
    {"s": "E. coli AIC222 (CRE)",           "g": "-", "m": True},
    {"s": "E. coli BAA-3170 (CRE)",         "g": "-", "m": True},
    {"s": "K. pneumoniae ATCC 13883",       "g": "-", "m": False},
    {"s": "K. pneumoniae BAA-2342 (ESBL)",  "g": "-", "m": True},
    {"s": "P. aeruginosa PAO1",             "g": "-", "m": False},
    {"s": "P. aeruginosa PA14",             "g": "-", "m": False},
    {"s": "P. aeruginosa BAA-3197 (FQR)",   "g": "-", "m": True},
    {"s": "S. enterica ATCC 9150",          "g": "-", "m": False},
    {"s": "S. enterica Typhimurium 700720", "g": "-", "m": False},
    {"s": "B. subtilis ATCC 23857",         "g": "+", "m": False},
    {"s": "S. aureus ATCC 12600",           "g": "+", "m": False},
    {"s": "S. aureus BAA-1556 (MRSA)",      "g": "+", "m": True},
    {"s": "E. faecalis 700802 (VRE)",       "g": "+", "m": True},
    {"s": "E. faecium 700221 (VRE)",        "g": "+", "m": True},
    {"s": "E. coli K-12 BW25113",           "g": "-", "m": False},
]

NUM_STRAINS = len(STRAIN_INFO)


# ══════════════════════════════════════════════════════════════════════════
# Data loading — all loaders strip whitespace from peptide names
# ══════════════════════════════════════════════════════════════════════════

def _name(row):
    return row['short_name'].strip()


def _mean(vals):
    return round(sum(vals) / len(vals), 1) if vals else None


def load_reference_table(path):
    with open(path) as f:
        return {_name(r): r for r in csv.DictReader(f)}


def load_mic(path):
    mic = {}
    with open(path) as f:
        rd = csv.DictReader(f)
        cols = [c for c in rd.fieldnames if c != 'short_name']
        for row in rd:
            mic[_name(row)] = [float(row[c].strip()) if row[c].strip() else None
                               for c in cols]
    return mic


def load_single_col(path, col):
    data = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                data[_name(row)] = min(float(row[col]), 1e6)
            except (ValueError, KeyError):
                pass
    return data


def load_membrane_assay(path):
    """Load DiSC3(5) or NPN: {name: [MaxRel, AUC]}."""
    data = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                data[_name(row)] = [float(row['MaxRel']), float(row['AUC'])]
            except (ValueError, KeyError):
                pass
    return data


def load_bestsel(path):
    """Returns {name: [[fH, fBa, fBp, fT, fO] × 4 solvents]}."""
    fracs = ['fH', 'f_beta_anti', 'f_beta_par', 'fturn', 'fothers']
    solvents = ['H2O', 'MeOH_H2O', 'SDS_H2O', 'TFE_H2O']
    data = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            data[_name(row)] = [
                [round(float(row.get(f'{fr}_{sv}', 0) or 0), 1) for fr in fracs]
                for sv in solvents
            ]
    return data


def load_lps(path):
    """Returns {name: [bc_2uM, bc_8uM, bc_32uM, bc_32uM_ca]}."""
    raw = defaultdict(lambda: defaultdict(list))
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                v = float(row['BC_displacement'])
            except (ValueError, KeyError):
                continue
            raw[_name(row)][(row['concentration'].strip(),
                             row['calcium'].strip())].append(v)
    keys = [('2', 'no'), ('8', 'no'), ('32', 'no'), ('32', 'yes')]
    return {sn: [_mean(grp.get(k, [])) for k in keys]
            for sn, grp in raw.items()}


def load_dna(path):
    """Returns {name: [deltaAmax, deltaA_12.5, deltaA_25, deltaA_50]}."""
    raw = defaultdict(dict)
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                v = float(row['value'])
            except (ValueError, KeyError):
                continue
            key = f"{row['metric'].strip()}@{row['concentration'].strip()}"
            raw[_name(row)][key] = v
    lookup = ['deltaAmax@', 'deltaA@12.5', 'deltaA@25.0', 'deltaA@50.0']
    return {sn: [round(v.get(k, 0), 2) for k in lookup]
            for sn, v in raw.items()}


def _load_if_exists(path, loader):
    return loader(path) if os.path.exists(path) else {}


def load_help(help_dir):
    """Load all .md files from help_dir into a dict keyed by stem name."""
    texts = {}
    if not os.path.isdir(help_dir):
        return texts
    for fname in sorted(os.listdir(help_dir)):
        if fname.endswith('.md'):
            key = fname[:-3]  # strip .md
            with open(os.path.join(help_dir, fname)) as f:
                texts[key] = f.read()
    return texts


# ══════════════════════════════════════════════════════════════════════════
# Build peptide list
# ══════════════════════════════════════════════════════════════════════════

def build_peptides(ref, mic_data, cc50_data, hc50_data,
                   disc_data, npn_data, bestsel_data, lps_data, dna_data):
    strain_mdr = [s['m'] for s in STRAIN_INFO]
    peptides = []

    for sn, r in ref.items():
        seq = r.get('sequence', '')
        mics = mic_data.get(sn, [None] * NUM_STRAINS)
        if len(mics) < NUM_STRAINS:
            mics += [None] * (NUM_STRAINS - len(mics))

        cc = cc50_data.get(sn)
        hc = hc50_data.get(sn)
        gmean = geometric_mean(mics)
        ns = sum(1 for v in mics if v is not None)
        nm = sum(1 for i, v in enumerate(mics) if v is not None and strain_mdr[i])

        pt = {
            'n': sn, 's': seq, 'c': r['category'],
            'm': r.get('conditioning', ''),
            'p': r.get('prototype_display') or r.get('prototype', ''),
            'o': r.get('objective', ''),
            'L': int(r['length']) if r.get('length') else 0,
            'mics': mics, 'mic': gmean,
            'cc': round(cc, 3) if cc else None,
            'hc': round(hc, 3) if hc else None,
            'ns': ns, 'nm': nm,
            'disc': disc_data.get(sn), 'npn': npn_data.get(sn),
            'bsl': bestsel_data.get(sn), 'lps': lps_data.get(sn),
            'dna': dna_data.get(sn),
        }
        pt.update(compute_descriptors(seq))
        pt['ti'] = round(cc / gmean, 2) if (cc and gmean) else None
        pt['th'] = round(hc / gmean, 2) if (hc and gmean) else None
        peptides.append(pt)

    return peptides


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description='Build OmegAMP dashboard')
    ap.add_argument('--data-dir', default='.', help='CSV directory (default: .)')
    ap.add_argument('--template', default=None, help='Template HTML path')
    ap.add_argument('--output', '-o', default='omegamp_dashboard.html')
    ap.add_argument('--help-dir', default=None,
                    help='Directory with .md help files (default: <data-dir>/help)')
    args = ap.parse_args()
    dd = args.data_dir

    tpl = args.template or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'dashboard_template.html')
    if not os.path.exists(tpl):
        sys.exit(f"ERROR: Template not found: {tpl}")

    print(f"Loading data from {dd}/")
    ref  = load_reference_table(os.path.join(dd, 'omegamp_reference_table.csv'))
    mic  = load_mic(os.path.join(dd, 'mic.csv'))
    cc50 = load_single_col(os.path.join(dd, 'cc50.csv'), 'CC50')
    hc50 = load_single_col(os.path.join(dd, 'hc50.csv'), 'HC50')
    disc = _load_if_exists(os.path.join(dd, 'disc.csv'), load_membrane_assay)
    npn  = _load_if_exists(os.path.join(dd, 'npn.csv'), load_membrane_assay)
    bsl  = _load_if_exists(os.path.join(dd, 'bestsel.csv'), load_bestsel)
    lps  = _load_if_exists(os.path.join(dd, 'lps_binding.csv'), load_lps)
    dna  = _load_if_exists(os.path.join(dd, 'dna_binding.csv'), load_dna)

    help_dir = args.help_dir or os.path.join(dd, 'help')
    help_texts = load_help(help_dir)
    if help_texts:
        print(f"  Loaded {len(help_texts)} help files: {', '.join(help_texts.keys())}")

    sizes = {'ref': len(ref), 'MIC': len(mic), 'CC50': len(cc50), 'HC50': len(hc50),
             'DiSC': len(disc), 'NPN': len(npn), 'BeStSel': len(bsl),
             'LPS': len(lps), 'DNA': len(dna)}
    print(f"  {' · '.join(f'{v} {k}' for k, v in sizes.items())}")

    print("  Computing physicochemical descriptors (modlAMP)...")
    peptides = build_peptides(ref, mic, cc50, hc50, disc, npn, bsl, lps, dna)
    cats = Counter(p['c'] for p in peptides)
    n_mic = sum(1 for p in peptides if any(v is not None for v in p['mics']))
    print(f"  Built {len(peptides)} peptides: {dict(cats)}")
    print(f"  {n_mic} with MIC data")

    payload = json.dumps({'peptides': peptides, 'strains': STRAIN_INFO,
                          'help': help_texts},
                         separators=(',', ':'))
    print(f"  JSON: {len(payload)/1024:.1f} KB")

    with open(tpl) as f:
        template = f.read()
    assert template.count('__DATA__') == 1, "Template must have exactly one __DATA__"
    html = template.replace('__DATA__', payload, 1)

    with open(args.output, 'w') as f:
        f.write(html)
    print(f"✓ {args.output} ({os.path.getsize(args.output)/1024:.0f} KB)")


if __name__ == '__main__':
    main()
