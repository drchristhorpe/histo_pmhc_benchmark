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
- `cutoff.py` — training-cutoff analysis: release-date loading (local CSV +
  RCSB/PDBe fallback), before/after classification, Mann-Whitney summaries, and
  the canonical pseudosequence-position index for per-residue figures.
- `pipeline.py` — `run_pipeline(...)` orchestrates all stages and writes outputs.
- `figures.py` — the 4-panel summary figure, the whole-structure cutoff box plots,
  and the per-method per-residue cutoff box plots.
- `cli.py` — Click CLI (`run`, `index`, `gaps`, `cutoff`).

## Training cutoffs

`Method.training_cutoff` holds each method's PDB cutoff. **HistoFold = 2018-04-30
(AlphaFold2), NOT AlphaFold3's 2021-09-30** — it builds on AlphaFold2. A test
enforces this. Per-residue interface figures aggregate by the canonical
pseudosequence rank but label the x-axis with the real MHC residue numbers
(`resnum` itself varies across structures, so it can't be the grouping key).

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
