---
name: histo-pmhc-benchmark
description: >
  Benchmark peptide-MHC class I structure predictions (AlphaFold3, Boltz-2,
  ESMFold2, ESMFold2-fast, HistoFold, or new methods) against ground-truth
  structures. Use when you have a folder of predictions and aligned ground-truth
  PDBs and need peptide RMSD (whole/per-residue/per-atom), peptide confidence,
  MHC binding-interface RMSD (mainchain + sidechain) and confidence at
  NetMHCpan-4.1 pseudosequence positions, per-method missing-complex gap CSVs,
  and a summary figure — as reproducible tables in one command.
---

# histo-pmhc-benchmark

Benchmarks pMHC class I structure predictions against ground truth. Composes
`histo_aligner` (common-frame superposition) and `histo_pseudosequence`
(NetMHCpan-4.1 interface positions). Ships as a Python library + `histo-pmhc-benchmark` CLI.

## When to use

The user has a dataset directory:

```
DATASET/
  predictions/<method>/<complex folders>/
  ground_truth/aligned/<pdb>_<copy>.pdb
  metadata/unique_structures.json   # optional
```

and wants the peptide + MHC-interface accuracy/confidence analysis, or just a
coverage/gap check. Also use it to **re-run** the analysis when a fuller set of
predictions lands — it globs the folders, so new complexes are picked up.

## CLI

- **Full analysis**: `histo-pmhc-benchmark run DATASET -o results/`
- **Re-score without re-aligning** (reuse aligned PDBs):
  `histo-pmhc-benchmark run DATASET -o results/ --aligned-root results/aligned --no-align`
- **Coverage table**: `histo-pmhc-benchmark index DATASET`
- **Missing-complex CSVs only**: `histo-pmhc-benchmark gaps DATASET -o results/`

## Library / kernel

```python
from histo_pmhc_benchmark import run_pipeline
manifest = run_pipeline(base="DATASET", out_dir="results", do_align=True, do_figure=True)
manifest["headline"]   # per-method n + median peptide Cα / interface RMSDs
```

Individual stages (`histo_pmhc_benchmark.stages`) run on a `build_dataset(base)`
`Dataset` for custom slices.

## Outputs

CSVs in `OUT/`: peptide RMSD (`_whole/_per_residue/_per_atom`), peptide confidence,
interface RMSD, interface confidence, two `*_calibration.csv`, `missing_<method>.csv`,
plus `pmhc_summary_figure.png` and `pmhc_run_manifest.json`.

## Interpreting

- Peptide **Cα RMSD** = backbone placement; **all-atom** adds rotamers. Per-residue
  shows the anchored-termini / bulged-centre pattern.
- Interface **mainchain** is near-solved (~0.2 Å); **sidechain** (0.5-1.2 Å) is the
  frontier that governs peptide specificity.
- Calibration **Spearman ρ** is negative when confidence tracks accuracy (higher
  pLDDT → lower RMSD); near zero = uncalibrated.

## Notes

- RMSD is direct-frame (no re-superposition) — relies on `histo_aligner` putting
  everything in the 1hhk frame first; the `run` command does this unless `--no-align`.
- Add a new prediction method by editing `config.py::METHODS` (one entry) — all
  stages, gaps, figure and colours follow.
- Method colours use Okabe-Ito, deliberately disjoint from the IBM bin/cutoff
  palette this project reserves for quality bins and training-cutoff semantics.
