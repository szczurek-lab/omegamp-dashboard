## DiSC3(5) -- Cytoplasmic Membrane Depolarization

The DiSC3(5) (3,3'-dipropylthiadicarbocyanine iodide) assay measures the integrity of the cytoplasmic (inner) membrane by monitoring membrane potential. DiSC3(5) accumulates in polarized membranes and self-quenches; when a peptide dissipates the membrane potential, the dye is released, producing a fluorescence increase.

### Gram-negative vs. Gram-positive context

The cytoplasmic membrane is the shared target across both Gram-negative and Gram-positive bacteria, but the path to reaching it differs. In Gram-negative bacteria, AMPs must first cross the LPS-containing outer membrane before encountering the cytoplasmic membrane. In Gram-positive bacteria, the thick peptidoglycan layer is porous to small cationic peptides, giving AMPs more direct access to the cytoplasmic membrane. This assay is measured against the Gram-negative indicator strain *A. baumannii*, so high DiSC3(5) activity implies that the peptide successfully traverses the outer membrane and then depolarizes the inner membrane. Peptides with high DiSC3(5) activity against Gram-negatives often also show activity against Gram-positive species, consistent with cytoplasmic membrane disruption being a shared vulnerability, though Gram-positive activity also correlates with greater hemolytic risk.

### Experimental conditions

All measurements were performed against *A. baumannii* ATCC 19606 at each peptide's MIC. Triton X-100 (complete membrane lysis) serves as the positive control; Polymyxin B and Levofloxacin are included as antibiotic references.

### Reported metrics

- **MaxRel (%)**: Maximum fluorescence intensity relative to the untreated baseline. Higher values indicate greater cytoplasmic membrane depolarization at peak response.
- **AUC**: Area under the fluorescence-time curve (integrated response). Captures both intensity and duration of depolarization.

### Reading the scatter plot

Peptides in the **upper-right quadrant** (high MaxRel, high AUC) are strong depolarizers with sustained cytoplasmic membrane disruption. **Upper-left** (low MaxRel, high AUC) are slow responders. **Lower-right** (high MaxRel, low AUC) are transient responders. Dashed lines mark the median across all tested peptides.

### Mechanistic context

DiSC-dominant peptides (high DiSC3(5), low NPN) disrupt the cytoplasmic membrane without strongly permeabilizing the outer membrane, suggesting they reach the inner membrane through mechanisms other than gross outer membrane lysis (e.g., self-promoted uptake or transient pore formation). DiSC-dominant activity is more strongly associated with hemolytic toxicity than NPN-dominant activity, likely because the cytoplasmic membrane disruption mechanism also affects eukaryotic cell membranes, which share structural similarity with bacterial cytoplasmic membranes.
