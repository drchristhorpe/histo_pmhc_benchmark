# histo-pmhc-benchmark

Benchmark **peptide-MHC class I structure predictions** against ground-truth
structures — peptide accuracy, MHC binding-groove accuracy, and confidence
calibration — for multiple prediction methods in one command.

It composes two histo_tools ([`histo_aligner`](https://github.com/drchristhorpe/histo_aligner)
for common-frame superposition and
[`histo_pseudosequence`](https://github.com/drchristhorpe/histo_pseudosequence)
for the NetMHCpan-4.1 interface positions) and ships as:

- a Python library — `from histo_pmhc_benchmark import run_pipeline`
- a CLI tool — `histo-pmhc-benchmark`
- a [Claude skill](skills/histo-pmhc-benchmark/SKILL.md)

Requires Python 3.14+.

## What it computes

For each prediction, aligned into the common 1hhk frame:

1. **Peptide RMSD** — whole peptide (Cα + all-atom), per-residue, per-atom, vs the
   best-matching ground-truth copy (direct-frame, no re-superposition).
2. **Peptide confidence** — per-residue / per-atom pLDDT (ESMFold b-factors
   rescaled to 0-100), plus a confidence↔accuracy calibration.
3. **MHC interface RMSD** — mainchain vs sidechain, at the 34 NetMHCpan-4.1
   pseudosequence positions.
4. **MHC interface confidence** — mainchain vs sidechain pLDDT at those positions,
   plus calibration.
5. **Gap analysis** — per-method CSVs of ground-truth complexes not predicted, in
   the `pdb_code,locus,allele_slug,peptide_sequence,resolution` schema.
6. **A 4-panel summary figure**.

## Install

```bash
uv sync                                  # dev environment from a checkout
uv sync --extra tools                    # also pull histo_aligner + histo_pseudosequence from GitHub
uv tool install .                        # install the CLI globally
```

The two composed histo_tools are declared as dependencies; in a dev checkout that
already has them on the path (e.g. `--system-site-packages`), they resolve from
there. `--extra tools` installs them from GitHub.

## CLI usage

```
histo-pmhc-benchmark run DATASET  [-o OUT] [--aligned-root DIR] [--no-align] [--no-figure] [--quiet]
histo-pmhc-benchmark index DATASET
histo-pmhc-benchmark gaps  DATASET [-o OUT]
```

`DATASET` is a directory containing:

```
DATASET/
  predictions/<method>/<complex folders>/   # alphafold3, boltz2, esmfold2-2026-05, esmfold2-fast-2026-05, histofold
  ground_truth/aligned/<pdb>_<copy>.pdb      # pre-aligned GT (1hhk frame)
  metadata/unique_structures.json            # optional allele/peptide backfill
```

- **Full run** (align + score + figure):
  `histo-pmhc-benchmark run DATASET -o results/`
- **Re-score without re-aligning** (reuse existing aligned PDBs):
  `histo-pmhc-benchmark run DATASET -o results/ --aligned-root results/aligned --no-align`
- **Quick coverage check**: `histo-pmhc-benchmark index DATASET`
- **Just the missing-complex CSVs**: `histo-pmhc-benchmark gaps DATASET -o results/`

Outputs (in `OUT/`): `pmhc_peptide_rmsd_{whole,per_residue,per_atom}.csv`,
`pmhc_peptide_confidence_*.csv`, `pmhc_interface_rmsd_*.csv`,
`pmhc_interface_confidence_*.csv`, `*_calibration.csv`, `missing_<method>.csv`,
`pmhc_summary_figure.png`, and `pmhc_run_manifest.json` (file map + headline
medians).

## Library usage

```python
from histo_pmhc_benchmark import run_pipeline

manifest = run_pipeline(
    base="path/to/DATASET",
    out_dir="results",
    do_align=True,        # False to reuse aligned_root
    do_figure=True,
)
manifest["headline"]      # {method: {n, peptide_ca_median, interface_mainchain_median, ...}}
manifest["files"]         # {logical_name: path}
```

Individual stages are in `histo_pmhc_benchmark.stages` and can be called on a
`Dataset` from `build_dataset(base)` for custom workflows.

## Adding a method

Prediction methods and their file/confidence adapters live in one place —
`config.py::METHODS`. Add a `Method(key, label, slug, plddt_scale)` entry and, if
its top-model path differs, extend `Method.top_model_path`. Everything downstream
(scoring, gaps, figure, colours) picks it up automatically.

## Development

```bash
uv sync --extra tools && uv run pytest
```

Tests exercise the geometry, config, gap-schema and calibration logic on small
inputs; the scoring stages are validated end-to-end against the real dataset. See
`docs/PLAN.md` for design rationale and `CHANGELOG.md` for release history.
