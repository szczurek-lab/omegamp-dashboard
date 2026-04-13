## Secondary Structure -- Circular Dichroism / BeStSel

Circular dichroism (CD) spectroscopy measures the differential absorption of left- and right-circularly polarized light by chiral molecules. Different secondary structure elements produce characteristic CD spectra, allowing estimation of the fractional composition of a peptide's structure.

### BeStSel analysis

Spectra are deconvoluted using BeStSel (Beta Structure Selection), an algorithm that estimates the fractional content of secondary structure elements from far-UV CD data. Five fractions are reported:

- **alpha-Helix** (red): The amphipathic alpha-helix is the dominant structural motif for membrane-active AMPs. Helix formation upon membrane contact segregates hydrophobic and cationic residues onto opposing faces, enabling membrane insertion and disruption.
- **beta-sheet antiparallel** (blue): Beta-sheet structures can indicate oligomerization or aggregation. Some AMPs form beta-sheet pores in membranes.
- **beta-sheet parallel** (light blue): Less common in short peptides; may indicate intermolecular association.
- **Turn** (amber): Beta-turns and other turn structures connecting secondary structure elements.
- **Other/disordered** (gray): Unstructured or irregular conformations.

### Solvent conditions

Structure is measured in four environments to probe conformational changes:

- **H2O**: Aqueous buffer (baseline). Most AMPs are unstructured in water.
- **MeOH/H2O**: Methanol-water mixture. Promotes helical folding by reducing solvent polarity; a mild membrane-mimetic environment.
- **SDS/H2O**: Sodium dodecyl sulfate micelles. The primary membrane-mimetic condition. Helicity in SDS indicates the peptide adopts an amphipathic helix upon membrane contact, the mechanistically relevant conformation for membranolytic AMPs.
- **TFE/H2O**: Trifluoroethanol-water mixture. A strong helix-inducing solvent that reveals the maximum helical propensity of the sequence.

### Reading the dashboard

The stacked bars show the fractional composition across all four solvents. The key comparison is between H2O (unstructured baseline) and SDS (membrane-mimetic). Peptides that transition from disordered in water to helical in SDS demonstrate the conformational switch characteristic of functional AMPs. The "% alpha" label shows the alpha-helix fraction for quick comparison.
