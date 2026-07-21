"""Direct-frame RMSD helpers.

All structures are pre-superposed onto the common 1hhk frame, so RMSDs are
computed atom-for-atom by residue-number correspondence with NO re-superposition.
"""

from __future__ import annotations

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

_PDB = PDBParser(QUIET=True)
_CIF = MMCIFParser(QUIET=True)


def load_model(path: str):
    """Parse the first model of a .pdb/.ent or .cif/.mmcif file."""
    p = path.lower()
    parser = _CIF if p.endswith((".cif", ".mmcif")) else _PDB
    return parser.get_structure("s", path)[0]


def chain_atoms(model, chain_id: str, heavy_only: bool = True) -> dict[int, dict[str, np.ndarray]]:
    """Map residue number -> {atom name: coord} for standard residues of a chain."""
    out: dict[int, dict[str, np.ndarray]] = {}
    if chain_id not in [c.id for c in model]:
        return out
    for r in model[chain_id]:
        if r.id[0] != " ":
            continue
        atoms = {a.name: a.coord for a in r if (not heavy_only or a.element != "H")}
        if atoms:
            out[r.id[1]] = atoms
    return out


def rmsd_over(gt: dict, pred: dict, positions=None, atomset=None, exclude=None):
    """Direct-frame RMSD over matched atoms.

    positions: residue numbers to include (None = all shared residues).
    atomset:  restrict to these atom names (e.g. mainchain).
    exclude:  drop these atom names (e.g. mainchain -> sidechain).
    Returns (rmsd, n_atoms); (nan, 0) when nothing matches.
    """
    keys = set(gt) & set(pred)
    if positions is not None:
        keys &= {p for p in positions if p is not None}
    sq = []
    for rn in keys:
        names = set(gt[rn]) & set(pred[rn])
        if atomset is not None:
            names &= set(atomset)
        if exclude is not None:
            names -= set(exclude)
        for an in names:
            d = gt[rn][an] - pred[rn][an]
            sq.append(float(np.dot(d, d)))
    if not sq:
        return float("nan"), 0
    return float(np.sqrt(np.mean(sq))), len(sq)


def ca_dev(gt_res: dict, pred_res: dict):
    """Cα deviation (Å) for one matched residue, or nan."""
    if "CA" in gt_res and "CA" in pred_res:
        return float(np.linalg.norm(gt_res["CA"] - pred_res["CA"]))
    return float("nan")
