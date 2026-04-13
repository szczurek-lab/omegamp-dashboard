# OmegAMP Interactive Dashboard

Interactive companion to the OmegAMP paper -- a generative diffusion framework for antimicrobial peptide design.

## Quick start

```bash
pip install modlamp

# Build locally
python build_dashboard.py --data-dir data --output docs/index.html

# Or build + encrypt with password
npm install -g staticrypt   # only needed for password protection
./deploy.sh --password YOUR_SECRET_HERE
```

## Deploy with GitHub Actions (recommended)

The repository includes a workflow that builds the dashboard and deploys it to GitHub Pages automatically on every push to `main`.

### Setup

1. Push this repo to GitHub.
2. Go to **Settings > Pages > Source** and select **GitHub Actions**.
3. (Optional) To enable password protection, go to **Settings > Secrets and variables > Actions** and add a repository secret named `DASHBOARD_PASSWORD`.
4. Push to `main` or trigger the workflow manually from the **Actions** tab.

The workflow will:
- Install Python + modlAMP
- Run `build_dashboard.py` to produce the HTML
- Encrypt with StatiCrypt if `DASHBOARD_PASSWORD` is set
- Deploy to GitHub Pages

If no password secret is configured the dashboard is published unprotected.

### Manual deploy (alternative)

For environments without GitHub Actions, the `deploy.sh` script builds the dashboard into `docs/` for branch-based Pages deployment:

```bash
./deploy.sh --password YOUR_SECRET_HERE
git add -A && git commit -m 'update dashboard' && git push
```

Then set **Settings > Pages > Source** to deploy from branch `main`, folder `/docs`.

## Password protection

The dashboard is AES-encrypted with [StatiCrypt](https://github.com/robinmoisson/staticrypt). Without the password, the page source is an encrypted blob -- peptide sequences, MIC data, and all assay results are unreadable.

- Visitors see a login page and must enter the password
- "Remember me" saves a hashed credential in the browser for 30 days
- To change the password: update the `DASHBOARD_PASSWORD` secret and re-run the workflow

## Views

| Tab | Description |
|-----|-------------|
| **Activity** | MIC across 20 strains + CC50/HC50/TI, sortable by family, group, or any column |
| **Physicochemical Landscape** | 8 descriptors (charge, hydrophobicity, moment, aromaticity, ...) with strain-selectable MIC coloring |
| **Membrane Activity** | NPN/DiSC3(5) quadrant plots for membrane disruption profiling |
| **Target Binding** | LPS neutralization (dose-response +/- Ca2+) and DNA binding (deltaA), family-grouped bar charts |
| **Structure** | BeStSel secondary structure across 4 solvents |
| **Safety** | MIC vs CC50/HC50 scatter with TI diagonal lines, strain selector |

## Features

- **Experiment group filters** -- De novo, Inactive-to-active, LPS-binding, DNA-binding with All/None buttons
- **Shared strain selector** -- choose strains for MIC summary across Safety and Physicochemical views
- **MIC aggregation toggle** -- switch between geometric mean, MIC50 (median), and MIC90 (90th percentile) globally across all views
- **Sub-1 μM MIC support** -- color scales and display extend to 0.125 μM for high-potency measurements
- **Contextual help panels** -- collapsible "?" buttons next to assay sections render markdown descriptions from `help/*.md` files
- Shared family/category filters across all tabs
- Search by peptide name
- Pin peptides -- persists across tab switches
- Brush-to-zoom on scatter plots
- Click column headers to sort in heatmap
- Export CSV -- download currently filtered data with selected MIC aggregation
- URL state -- shareable links that restore exact view

## Input data

| File | Rows | Required | Description |
|------|------|----------|-------------|
| `omegamp_reference_table.csv` | 217 | yes | Peptide metadata |
| `mic.csv` | 217 | yes | MIC across 20 bacterial strains (supports sub-1 μM values) |
| `cc50.csv` | 215 | yes | Cytotoxicity CC50 |
| `hc50.csv` | 215 | yes | Hemolysis HC50 |
| `disc.csv` | 196 | | DiSC3(5) membrane depolarization |
| `npn.csv` | 196 | | NPN outer membrane permeabilization |
| `bestsel.csv` | 190 | | BeStSel secondary structure |
| `lps_binding.csv` | 577 | | LPS binding (BC displacement) |
| `dna_binding.csv` | 111 | | DNA binding (deltaA) |
| `help/*.md` | -- | | Markdown help files (see below) |

Optional files are loaded if present; tabs with missing data render empty.

### Help files

Place `.md` files in a `help/` subdirectory (or pass `--help-dir`). The filename stem must match the `data-help` attribute of a help button in the template:

| File | Section |
|------|---------|
| `lps-binding.md` | LPS neutralization (Target Binding tab) |
| `dna-binding.md` | DNA binding (Target Binding tab) |
| `npn.md` | NPN assay (Membrane Activity tab) |
| `disc.md` | DiSC3(5) assay (Membrane Activity tab) |
| `bestsel.md` | Secondary structure (Structure tab) |

If a file is missing, its "?" button is hidden automatically.

## Strain panel

20 strains (15 Gram-, 5 Gram+), including 8 MDR clinical isolates:

| Resistance | Strain | ATCC # |
|------------|--------|--------|
| CRAB | *A. baumannii* | BAA-1605 |
| CRE | *E. coli* | AIC222, BAA-3170 |
| ESBL | *K. pneumoniae* | BAA-2342 |
| FQR | *P. aeruginosa* | BAA-3197 |
| MRSA | *S. aureus* | BAA-1556 |
| VRE | *E. faecalis*, *E. faecium* | 700802, 700221 |

## Dependencies

- **Python 3.8+** with [modlAMP](https://modlamp.org/) -- physicochemical descriptors
- **Node.js** with [StatiCrypt](https://github.com/robinmoisson/staticrypt) -- password protection (optional)
- **D3.js v7** -- loaded from CDN at runtime, no local install needed
