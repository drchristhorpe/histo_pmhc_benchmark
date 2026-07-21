# CLAUDE.md — histo_pmhc_benchmark

## What this is

A pipeline that benchmarks pMHC class I structure predictions (multiple methods)
against ground-truth structures. Composes `histo_aligner` (common-frame
superposition onto 1hhk) and `histo_pseudosequence` (NetMHCpan-4.1 interface
positions). Ships as library + Click CLI + skill.

## Module map

- `config.py` — `METHODS` (the single source of truth for method definitions and
  per-method file/confidence adapters), colours, chain/atom constants. **Add a
  method here, nowhere else.**
- `dataset.py` — `build_dataset(base)` scans a dataset dir into a `Dataset`
  (GT copies, per-method complex folders, allele/peptide metadata).
- `geometry.py` — direct-frame RMSD helpers. Structures are PRE-ALIGNED into a
  common frame, so RMSD is atom-for-atom by residue number with **no
  re-superposition**. Do not add a superposition step.
- `stages.py` — the independently-callable stages (align, pseudosequence,
  peptide RMSD/confidence, interface RMSD/confidence, gaps, calibration).
- `pipeline.py` — `run_pipeline(...)` orchestrates all stages and writes outputs.
- `figures.py` — the 4-panel summary figure.
- `cli.py` — Click CLI (`run`, `index`, `gaps`).

## Key invariants

- **Direct-frame RMSD.** Everything relies on GT and predictions sharing the 1hhk
  frame (MHC 1-180 COM invariant). `histo_aligner` provides that; never re-fit.
- **Best-GT-copy selection.** Many complexes have multiple GT copies; each metric
  compares against the lowest-RMSD copy. Keep that convention.
- **ESMFold pLDDT scale.** ESMFold PDBs store pLDDT on 0-1 in the b-factor column;
  every other method uses 0-100. `Method.plddt_scale` handles this — don't
  hardcode a scale in the stages.
- **Method colours must NOT reuse the reserved IBM bin/cutoff hexes.** This project
  reserves the IBM colour-blind-safe set (#785EF0/#DC267F/#FE6100/#FFB000/#648FFF)
  for confidence-quality bins and training-cutoff/overlay semantics. Method
  identity uses the Okabe-Ito palette instead. A test enforces no overlap.
- **Gap CSV schema is fixed**: `pdb_code,locus,allele_slug,peptide_sequence,resolution`
  — pdb_code lower, peptide upper, resolution 0.0. A test enforces it.

## Re-running on new predictions

The stages glob the prediction folders, so dropping new complexes into
`predictions/<method>/` and re-running picks them up. Use `--no-align` +
`--aligned-root` to re-score without re-aligning unchanged structures.

## Testing

`uv run pytest`. Unit tests cover geometry/config/gap/calibration on small
inputs. Full-dataset validation is a manual `run` against the real dataset (the
headline medians are the regression check).

## Scope

Benchmarking only — it measures and reports. It does not predict structures, run
NetMHCpan, or fetch data. Keep new options minimal and check before broadening.
