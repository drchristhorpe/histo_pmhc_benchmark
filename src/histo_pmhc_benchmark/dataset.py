"""Index a pMHC prediction dataset: ground-truth copies, per-method complex
folders, and allele/peptide metadata."""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field

from histo_pmhc_benchmark.config import METHODS

_PDB_RE = re.compile(r"^([0-9][a-z0-9]{3})", re.I)


def _pid_from_name(name: str) -> str | None:
    n = name[5:] if name.startswith("fold_") else name
    m = _PDB_RE.match(n)
    return m.group(1).lower() if m else None


def _parse_folder(name: str) -> tuple[str | None, str | None, str | None]:
    """Return (pdb_id, allele_slug, peptide) parsed from a prediction folder name."""
    n = name[5:] if name.startswith("fold_") else name
    if "__" in n:
        p = n.split("__")
        return (p[0].lower(),
                p[1] if len(p) > 1 else None,
                p[2].upper() if len(p) > 2 else None)
    m = re.match(r"^([0-9][a-z0-9]{3})_(.+)$", n, re.I)
    if not m:
        pid = _pid_from_name(name)
        return pid, None, None
    pid = m.group(1).lower()
    toks = m.group(2).split("_")
    return (pid,
            "_".join(toks[:-1]) if len(toks) > 1 else None,
            toks[-1].upper() if toks else None)


@dataclass
class Dataset:
    """A resolved view of one pMHC benchmark dataset directory."""

    base: str
    gt_files: dict[str, list[tuple[int, str]]]        # pid -> [(copy_no, filename), ...]
    method_folders: dict[str, dict[str, str]]          # method.key -> {pid: folder}
    pdb_meta: dict[str, tuple[str, str]]               # pid -> (allele_slug, peptide)

    @property
    def gt_ids(self) -> set[str]:
        return set(self.gt_files)

    def gt_aligned_dir(self) -> str:
        return os.path.join(self.base, "ground_truth", "aligned")

    def to_index(self) -> dict:
        return {"gt_files": {k: [list(t) for t in v] for k, v in self.gt_files.items()},
                "method_folders": self.method_folders}


def build_dataset(base: str) -> Dataset:
    """Scan a dataset directory and resolve GT copies, method folders, metadata."""
    base = os.path.abspath(os.path.expanduser(base))
    gt_dir = os.path.join(base, "ground_truth", "aligned")
    if not os.path.isdir(gt_dir):
        raise FileNotFoundError(f"ground_truth/aligned not found under {base}")

    gt_files: dict[str, list[tuple[int, str]]] = {}
    for f in os.listdir(gt_dir):
        m = re.match(r"^([0-9][a-z0-9]{3})_(\d+)\.pdb$", f)
        if m:
            gt_files.setdefault(m.group(1).lower(), []).append((int(m.group(2)), f))
    for pid in gt_files:
        gt_files[pid].sort()

    method_folders: dict[str, dict[str, str]] = {}
    pdb_meta: dict[str, tuple[str, str]] = {}
    for meth in METHODS:
        mdir = os.path.join(base, "predictions", meth.key)
        d: dict[str, str] = {}
        if os.path.isdir(mdir):
            for fld in os.listdir(mdir):
                if fld == ".DS_Store":
                    continue
                pid, al, pep = _parse_folder(fld)
                if pid:
                    d[pid] = fld
                    if al and pep and pid not in pdb_meta:
                        pdb_meta[pid] = (al.lower(), pep.upper())
        method_folders[meth.key] = d

    # backfill allele/peptide from metadata/unique_structures.json when present
    uniq_path = os.path.join(base, "metadata", "unique_structures.json")
    if os.path.exists(uniq_path):
        uniq = json.load(open(uniq_path))
        for key, pdbs in uniq.items():
            parts = key.rsplit("_", 1)
            al = parts[0]
            pep = parts[1] if len(parts) > 1 else None
            for p in pdbs:
                pl = p.lower()
                if pl not in pdb_meta and al and pep:
                    pdb_meta[pl] = (al.lower(), pep.upper())

    return Dataset(base=base, gt_files=gt_files, method_folders=method_folders, pdb_meta=pdb_meta)
