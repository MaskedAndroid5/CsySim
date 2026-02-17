CsySim 🦈🧬  — First dedicated mechanistic simulator for the Type I-F CRISPR-Csy complex (Pseudomonas aeruginosa). Models strict 5'-GG PAM, seed hypersensitivity (positions 1–8), segmental R-loop unzipping along the Csy3 backbone, ΔG-based biophysical scoring, and stochastic Monte Carlo propagation.


The world's first dedicated mechanistic simulator for the Type I-F CRISPR-Csy surveillance complex  
Pseudomonas aeruginosa — the classic seahorse-shaped Type I-F system.

While Cas9 receives the vast majority of attention and tooling, **Type I-F** is one of the most abundant CRISPR-Cas systems in nature — yet it has almost no dedicated computational simulators. CsySim changes that.

This project provides a from-scratch, literature-grounded simulator focused exclusively on the Csy complex (Csy1, Csy2, Csy3₆, Csy4 + crRNA), with emphasis on:

- Strict 5'-GG PAM recognition on the non-target strand  
- Hypersensitive PAM-proximal seed region (positions 1–8 must be perfect)  
- **Segmental, directional R-loop unzipping** along the six Csy3 subunits  
- Stochastic propagation modeled via Monte Carlo trials  
- Biophysical ΔG scoring (nearest-neighbor RNA:DNA hybrid stability) blended with mismatch-based Gaussian penalties  
- Structural mapping of mismatches to approximate Csy3 residue contacts  

Features

- v0.1 — Fast deterministic binding predictor (PAM + seed strictness + Gaussian distal penalty)  
- v0.2 — Full segmental Monte Carlo simulator (directional unzipping, stall probabilities, interference proxy)  
- v0.2.2 — CsyZipper
  - Nearest-neighbor ΔG scoring for RNA:DNA hybrids (initiation + terminal penalties)  
  - Blended ΔG/Gaussian per-segment probabilities  
  - Structural notes linking mismatches to approximate Csy3 thumb-loop / contact residues  
  - Batch CSV processing for screening many targets  
  - Reproducible runs via `--seed`  

--- Installation

```bash
git clone https://github.com/YOUR_USERNAME/CsySim.git
cd CsySim
# No external dependencies required (pure stdlib Python)
