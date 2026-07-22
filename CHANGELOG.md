# Changelog

## 0.2.0

Training-cutoff analysis.

- New `cutoff` stage and CLI subcommand: split each method's predictions into
  structures released before vs after that method's PDB training-data cutoff, and
  compare accuracy (RMSD) and confidence (pLDDT) across the boundary with a
  Mann-Whitney U test. Release dates come from the dataset's local metadata CSV
  with an RCSB/PDBe fallback (`--offline` to disable).
- `Method` gains a `training_cutoff` field. **HistoFold uses the AlphaFold2
  cutoff (2018-04-30)** — it builds on AlphaFold2, not AlphaFold3.
- New figures: `pmhc_confidence_cutoff_boxplots.png` (whole-structure pLDDT
  box plots before/after) and one `pmhc_perres_cutoff_<method>.png` per method
  with an after-cutoff set — per-residue box plots of peptide pLDDT + Cα deviation
  and binding-site sidechain pLDDT + RMSD, the binding-site axis labelled with the
  canonical NetMHCpan-4.1 MHC residue numbers.
- `run` gains `--no-cutoff` and `--offline`; new `cutoff` subcommand runs the
  analysis on already-scored tables.

## 0.1.0

Initial release.

- `run_pipeline` + Click CLI (`run`, `index`, `gaps`) benchmarking pMHC class I
  structure predictions against ground truth.
- Peptide RMSD (whole / per-residue / per-atom) and confidence, MHC interface
  RMSD and confidence (mainchain vs sidechain) at NetMHCpan-4.1 pseudosequence
  positions, confidence↔accuracy calibration, per-method gap CSVs, and a 4-panel
  summary figure.
- Composes `histo_aligner` and `histo_pseudosequence`.
- Five methods supported out of the box: AlphaFold3, Boltz-2, ESMFold2,
  ESMFold2-fast, HistoFold. New methods add via one `config.py::METHODS` entry.
- Method colours use the Okabe-Ito palette, kept disjoint from the IBM
  bin/cutoff hexes reserved elsewhere in the project.
