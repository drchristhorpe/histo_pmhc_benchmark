# histo-pmhc-benchmark — design notes

## Purpose

Turn the ad-hoc pMHC benchmark analysis into a reproducible, re-runnable pipeline:
given a folder of predictions from several methods and a folder of aligned
ground-truth structures, produce every peptide/interface accuracy and confidence
table plus a summary figure, in one command.

## Why direct-frame RMSD

`histo_aligner` superposes every structure onto the canonical 1hhk MHC frame
using the conserved α1/α2 framework. Because GT and predictions then share one
coordinate frame, RMSD is computed atom-for-atom by residue-number
correspondence with **no per-pair re-superposition**. This measures true
placement error (including rigid-body drift of the peptide/groove relative to the
MHC platform), which a per-pair Kabsch fit would hide. It also makes per-residue
and per-atom deviations directly comparable across methods.

## Interface = NetMHCpan-4.1 pseudosequence

The MHC binding-groove analysis is restricted to the 34 NetMHCpan-4.1
peptide-contacting positions, located per-structure by `histo_pseudosequence`
(alignment-based transfer, robust to numbering differences). Mainchain
(N,Cα,C,O) and sidechain RMSD are reported separately because backbone geometry
is near-solved while rotamer placement is the frontier — collapsing them hides
the part that governs specificity.

## Best-GT-copy convention

Many complexes are crystallized with several copies in the asymmetric unit. Each
prediction is scored against the copy that minimizes the relevant whole-region
RMSD, so a method is not penalized for matching one valid copy over another.

## Method adapters in one place

Methods differ in (a) where the top-ranked model file sits and (b) the pLDDT
scale in the b-factor column (ESMFold 0-1, others 0-100). Both are captured in
`config.py::Method`, so the scoring stages stay method-agnostic and adding a
method is a one-line change.

## Colour discipline

The parent project reserves the IBM colour-blind-safe palette for confidence
-quality bins and training-cutoff/overlay semantics. Method identity is a
different categorical dimension, so it uses the Okabe-Ito qualitative palette to
avoid implying a bin/cutoff meaning. A test asserts the two sets don't overlap.

## Stage independence

Each stage is a pure function on a `Dataset` returning DataFrames, so the CLI,
the kernel helper, and any custom analysis can call them individually.
`run_pipeline` is just the standard composition and the file-writing/manifest
layer.

## Re-runnability

Stages glob the prediction directories, so a later, fuller set of predictions is
handled by re-running — no code change. `--no-align` reuses existing aligned PDBs
so only scoring re-runs when structures are unchanged.

## Scope

Measurement and reporting only. Prediction, NetMHCpan execution, and data
fetching are out of scope by design.
