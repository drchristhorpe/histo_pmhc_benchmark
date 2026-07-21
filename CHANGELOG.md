# Changelog

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
