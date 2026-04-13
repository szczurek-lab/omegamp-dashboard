## DNA Binding -- CD Spectrometry

This assay measures the ability of peptides to bind double-stranded DNA using circular dichroism (CD) spectrometry. When a peptide binds DNA, conformational changes in the DNA-peptide complex alter the CD spectrum. The change in absorbance (deltaA) at characteristic wavelengths quantifies the extent of binding.

### Reported metric

- **deltaA**: Change in absorbance upon peptide-DNA interaction, measured by CD. Higher deltaA values indicate stronger DNA binding. A deltaA of 0 means no detectable binding; values exceeding the template (prototype) sequence indicate enhanced binding relative to the parent scaffold.

### Concentrations

Binding is assessed at three peptide concentrations (12.5, 25, and 50 uM) to capture concentration-dependent behavior. **deltaA_max** reports the highest deltaA observed across all tested concentrations. The dashboard primarily displays the 12.5 uM condition, where the comparison between analogs and their prototype is cleanest.

### Design context

The DNA-binding analogs are derived from the GCN4 bZIP domain (PDB: 1GD2), a leucine zipper transcription factor that binds DNA through its basic region. OmegAMP's motif-guided generation fixed the leucine/valine residues at the coiled-coil interface (heptad a and d positions) while optimizing surrounding residues for antimicrobial activity. The goal is simultaneous DNA binding and antimicrobial function.

### Mechanistic note

Peptides that combine high NPN with low DiSC3(5) may reach intracellular DNA targets via membrane translocation rather than lysis, consistent with a non-membranolytic antimicrobial mechanism.
