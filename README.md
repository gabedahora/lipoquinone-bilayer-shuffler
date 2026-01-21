# Description
Randomize lipoquinone (UQ8/UQOL8) positions in lipid bilayers for GROMACS replicas using MDAnalysis: keeps molecules intact across PBC, avoids water, and enforces steric/protein-distance constraints.

lipoquinone-bilayer-shuffler is a command-line Python tool to generate replica-ready membrane systems by randomizing UQ8 (I'll expand later to UQOL8) positions inside a lipid bilayer. The tool preserves .gro atom order for GROMACS compatibility, wraps residues as whole molecules to prevent PBC splitting, biases the quinone/quinol headgroup to the bilayer interface, and applies configurable distance constraints to avoid overlaps and prevent insertion into tight protein grooves while still permitting annular membrane occupancy.


# What does it do?
Randomize **lipoquinone** positions (currently **UQ8**, extendable to **UQOL8**) inside a lipid bilayer for **GROMACS replica generation**, using **MDAnalysis**.
This tool repositions UQ molecules while:

- Preserving **.gro atom ordering** (safe for GROMACS topology mapping)
- Avoiding **PBC splitting** by wrapping residues as a whole (molecules stay intact)
- Keeping lipoquinones out of **bulk water** using phosphate-plane Z constraints
- Avoiding **overlaps** with protein, lipids, and previously placed quinones
- Preventing insertion into tight **protein grooves**, while still allowing annular occupancy
- Providing reproducibility via `--seed`

---

## Features
- **Rigid-body repositioning** of each lipoquinone:
  - Random XY placement
  - Leaflet-aware Z placement (head near interface region)
  - Head→tail axis aligned toward the bilayer midplane
  - Random in-plane (azimuthal) rotation for diversity
- **Two-tier steric thresholds** (separate constraints):
  - Strict distance to **protein** and previously placed **UQ**
  - Softer distance to **lipids** to improve insertion success in dense membranes
- **Membrane confinement**:
  - Ensures all lipoquinone **heavy atoms** remain between phosphate planes, with a configurable keepout.


---

## Basic Usage
- **General usage**:
    _python lipoquinone_bilayer_shuffler.py -i step5_input.gro -o shuffled.gro --seed 1234_
    
- **One I used**: 
    _python lipoquinone_bilayer_shuffler.py -i step5_input.gro -o shuffled.gro --seed 1001 --min-clash-prot-nm 0.20 --min-clash-lip-nm 0.10 --ring-protein-excl-nm 0.80 --z-margin-nm 0.50 --z-keepout-nm 0.10 --z-head-band-nm 0.80 --max-attempts 50000_

---

## Requirements
- Python 3.9+ (tested with Python 3.12)
- `numpy`
- `MDAnalysis`

Install with:

```bash
pip install numpy MDAnalysis
