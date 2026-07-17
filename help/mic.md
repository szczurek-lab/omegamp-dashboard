## MIC -- Minimum Inhibitory Concentration

MIC is the lowest peptide concentration that fully prevents visible bacterial growth. It is the primary measure of antimicrobial potency: lower MIC means a more potent peptide. Each peptide is tested against a panel of 20 bacterial strains spanning 14 Gram-negative, 5 Gram-positive, and 1 laboratory reference (*E. coli* K-12 BW25113), including **8 multidrug-resistant (MDR) clinical isolates**.

### Assay method

MIC values were determined by broth microdilution in untreated 96-well microplates. Peptides were prepared as twofold serial dilutions and mixed 1:1 with LB medium containing bacteria at mid-logarithmic phase. The MIC was read as the lowest concentration that fully prevented visible growth after 24 h at 37 C. Each assay was performed independently in triplicate.

### Value range and censoring

Most values fall on the twofold dilution ladder (1, 2, 4, ... 64 uM). Peptides that showed no inhibition at the highest concentration tested are right-censored as **>64 uM** (shown as empty cells in the heatmap) and count as tested-but-inactive. Highly potent peptides were re-measured below 1 uM, so values such as 0.5, 0.3, and 0.25 uM also appear.

### Summarizing across strains

The header toggle controls how a peptide's 20 per-strain MICs are collapsed to a single number:

- **MIC50** (median, default): the middle value across the panel -- a robust summary of overall potency.
- **Geometric mean**: the geometric average across strains.
- **MIC90**: the concentration inhibiting 90% of strains -- emphasizes the harder-to-reach tail.

The "Active if MIC <=" dropdown sets the threshold that counts a strain as reached (default 4 uM); the per-family "strains active" counts and the sorting presets use this cutoff.

### MDR panel

The 8 MDR isolates are carbapenem-resistant *A. baumannii* (BAA-1605), carbapenem-resistant *E. coli* (AIC222, BAA-3170), ESBL-producing *K. pneumoniae* (BAA-2342), fluoroquinolone-resistant *P. aeruginosa* (BAA-3197), MRSA *S. aureus* (BAA-1556), and vancomycin-resistant *E. faecalis* (700802) and *E. faecium* (700221). Reaching these at MIC <= 2 uM is the strictest coverage bar and drove lead selection.
