# OmegAMP Interactive Dashboard

Interactive companion to the OmegAMP paper — a generative diffusion framework for antimicrobial peptide design.

## Quick start

```bash
pip install modlamp
npm install -g staticrypt   # only needed for password protection

# Build + encrypt
./deploy.sh --password YOUR_SECRET_HERE

# Or build without password
./deploy.sh
```

## Password protection

The dashboard is AES-encrypted with [StatiCrypt](https://github.com/robinmoisson/staticrypt). Without the password, the page source is an encrypted blob — peptide sequences, MIC data, and all assay results are unreadable.

- Visitors see a login page and must enter the password
- "Remember me" saves a hashed credential in the browser for 30 days
- To change the password: `./deploy.sh --password NEW_PASSWORD` and push

For local development, an unencrypted copy is kept at `docs/index_unprotected.html`.

## Deploy on GitHub Pages

1. Push this repo to GitHub
2. **Settings → Pages → Source**: Deploy from branch `main`, folder `/docs`
3. Share the URL + password with collaborators

## Views

| Tab | Description |
|-----|-------------|
| **Therapeutic Index** | MIC vs CC50/HC50 scatter with TI diagonal lines, strain selector |
| **Physicochemical** | 8 descriptors (charge, hydrophobicity, moment, aromaticity, ...) |
| **Activity Heatmap** | MIC across 20 strains + CC50/HC50/TI, sortable by any column |
| **Membrane** | NPN/DiSC3(5) quadrant plots for membrane disruption profiling |
| **Binding** | LPS neutralization (dose-response ± Ca²⁺) and DNA binding (ΔA) |
| **Structure** | BeStSel secondary structure across 4 solvents |

## Features

- Shared family/category filters across all tabs
- Search by peptide name
- Pin peptides — persists across tab switches
- Brush-to-zoom on scatter plots
- Click column headers to sort in heatmap
- ⬇ Export CSV — download currently filtered data
- URL state — shareable links that restore exact view (`#hm/fam=BoCo1/pin=Ω-AP-BoCo1-3`)

## Input data

| File | Rows | Required | Description |
|------|------|----------|-------------|
| `omegamp_reference_table.csv` | 217 | ✓ | Peptide metadata |
| `mic.csv` | 217 | ✓ | MIC across 20 bacterial strains |
| `cc50.csv` | 215 | ✓ | Cytotoxicity CC50 |
| `hc50.csv` | 215 | ✓ | Hemolysis HC50 |
| `disc.csv` | 196 | | DiSC3(5) membrane depolarization |
| `npn.csv` | 196 | | NPN outer membrane permeabilization |
| `bestsel.csv` | 190 | | BeStSel secondary structure |
| `lps_binding.csv` | 577 | | LPS binding (BC displacement) |
| `dna_binding.csv` | 111 | | DNA binding (ΔA) |

Optional files are loaded if present; tabs with missing data render empty.

## Strain panel

20 strains (15 Gram−, 5 Gram+), including 8 MDR clinical isolates:

| Resistance | Strain | ATCC # |
|------------|--------|--------|
| CRAB | *A. baumannii* | BAA-1605 |
| CRE | *E. coli* | AIC222, BAA-3170 |
| ESBL | *K. pneumoniae* | BAA-2342 |
| FQR | *P. aeruginosa* | BAA-3197 |
| MRSA | *S. aureus* | BAA-1556 |
| VRE | *E. faecalis*, *E. faecium* | 700802, 700221 |

## Dependencies

- **Python 3.8+** with [modlAMP](https://modlamp.org/) — physicochemical descriptors
- **Node.js** with [StatiCrypt](https://github.com/robinmoisson/staticrypt) — password protection (optional)
- **D3.js v7** — loaded from CDN at runtime, no local install needed
