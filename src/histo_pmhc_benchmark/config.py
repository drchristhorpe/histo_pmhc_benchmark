"""Method definitions and per-method file/confidence adapters.

Each prediction method stores its top-ranked model in a different place and its
per-residue confidence (pLDDT) on a different scale. This module centralizes
those differences so the pipeline stages stay method-agnostic.
"""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass


@dataclass(frozen=True)
class Method:
    """One prediction method: folder name, display label, and how to find its
    top-ranked model file and normalize its confidence."""

    key: str            # prediction subfolder name
    label: str          # display label used in all output tables
    slug: str           # filesystem-safe short id
    plddt_scale: float  # multiply b-factor by this to get pLDDT on 0-100
    training_cutoff: str  # PDB training-data cutoff (ISO date); structures released
    #                       on/after this were unseen at training time.

    def top_model_path(self, base: str, folder: str) -> str | None:
        """Absolute path to this method's top-ranked model for one complex,
        or None if absent. `folder` is the complex's prediction subdirectory."""
        d = os.path.join(base, "predictions", self.key, folder)
        if self.key == "alphafold3":
            g = glob.glob(f"{d}/*_model_0.cif")
            return g[0] if g else None
        if self.key == "boltz2":
            f = f"{d}/outputs/files/prediction/sample_0_predicted_structure.cif"
            return f if os.path.exists(f) else None
        if self.key in ("esmfold2-2026-05", "esmfold2-fast-2026-05"):
            g = glob.glob(f"{d}/*rank_1*prediction_1.pdb")
            return g[0] if g else None
        if self.key == "histofold":
            g = glob.glob(f"{d}/*_relaxed_rank_001_*.pdb")
            return g[0] if g else None
        return None


# The five benchmarked methods. ESMFold PDBs carry pLDDT on a 0-1 scale in the
# b-factor column; every other method uses 0-100, so ESMFold is scaled x100.
#
# training_cutoff = the PDB release-date boundary of each method's training data:
#   AlphaFold3            2021-09-30  (Abramson et al. 2024)
#   Boltz-2              2023-06-01  (Boltz-2 preprint)
#   ESMFold2 / -fast     2021-09-30  (ESM training snapshot)
#   HistoFold            2018-04-30  (built on AlphaFold2 — Jumper et al. 2021 —
#                                     whose PDB training set ends 30 Apr 2018;
#                                     NOT AlphaFold3's 2021 cutoff)
METHODS: tuple[Method, ...] = (
    Method("alphafold3", "AlphaFold3", "alphafold3", 1.0, "2021-09-30"),
    Method("boltz2", "Boltz-2", "boltz2", 1.0, "2023-06-01"),
    Method("esmfold2-2026-05", "ESMFold2", "esmfold2", 100.0, "2021-09-30"),
    Method("esmfold2-fast-2026-05", "ESMFold2-fast", "esmfold2_fast", 100.0, "2021-09-30"),
    Method("histofold", "HistoFold", "histofold", 1.0, "2018-04-30"),
)

METHOD_ORDER: tuple[str, ...] = tuple(m.label for m in METHODS)

# Method label -> distinct CVD-safe plotting colour. Uses the Okabe-Ito
# qualitative palette (method identity) deliberately kept SEPARATE from the
# IBM bin/cutoff hexes this project reserves for confidence-quality bins and
# training-cutoff/overlay semantics, so method colour never collides with those.
METHOD_COLORS: dict[str, str] = {
    "AlphaFold3": "#0072B2",     # blue
    "Boltz-2": "#009E73",        # bluish green
    "ESMFold2": "#E69F00",       # orange
    "ESMFold2-fast": "#56B4E9",  # sky blue
    "HistoFold": "#CC79A7",      # reddish purple
}

# Peptide chain, MHC heavy chain, and the mainchain atom set.
PEPTIDE_CHAIN = "C"
MHC_CHAIN = "A"
MAINCHAIN_ATOMS = frozenset({"N", "CA", "C", "O"})


def method_by_key(key: str) -> Method | None:
    for m in METHODS:
        if m.key == key:
            return m
    return None


def method_by_label(label: str) -> Method | None:
    for m in METHODS:
        if m.label == label:
            return m
    return None
