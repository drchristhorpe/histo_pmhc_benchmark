"""Pipeline stages. Each returns a pandas DataFrame (or dict of them) and is
independently callable; the CLI and kernel helper compose them."""

from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd

from histo_pmhc_benchmark.config import (
    MAINCHAIN_ATOMS,
    METHODS,
    METHOD_ORDER,
    MHC_CHAIN,
    PEPTIDE_CHAIN,
    method_by_label,
)
from histo_pmhc_benchmark.dataset import Dataset
from histo_pmhc_benchmark.geometry import ca_dev, chain_atoms, load_model, rmsd_over

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------- alignment ---
def align_predictions(ds: Dataset, out_root: str, progress=None) -> pd.DataFrame:
    """Superpose each method's top model onto the 1hhk frame via histo_aligner.

    Writes `<out_root>/<method.key>/<pid>_aligned.pdb`. Returns a per-structure
    report (status, align_rmsd, framework_matched, mhc_com).
    """
    from histo_aligner.core import HistoAligner, save_structure

    aln = HistoAligner()
    rows = []
    for meth in METHODS:
        folders = ds.method_folders.get(meth.key, {})
        outdir = os.path.join(out_root, meth.key)
        os.makedirs(outdir, exist_ok=True)
        n_ok = 0
        for pid, fld in folders.items():
            src = meth.top_model_path(ds.base, fld)
            if not src:
                rows.append({"method": meth.label, "pdb_id": pid, "status": "no_model"})
                continue
            try:
                res = aln.align(src)
                outpath = os.path.join(outdir, f"{pid}_aligned.pdb")
                save_structure(res.structure, outpath, format="pdb")
                mcoords = [a.coord for a in res.structure[0][MHC_CHAIN].get_atoms()
                           if a.get_parent().id[0] == " " and 1 <= a.get_parent().id[1] <= 180]
                com = list(np.round(np.mean(mcoords, axis=0), 3)) if mcoords else None
                rows.append({"method": meth.label, "pdb_id": pid, "status": "ok",
                             "align_rmsd": round(res.rmsd, 4),
                             "framework_matched": res.framework_matched_count,
                             "mhc_com_x": com[0] if com else None,
                             "mhc_com_y": com[1] if com else None,
                             "mhc_com_z": com[2] if com else None})
                n_ok += 1
            except Exception as e:  # noqa: BLE001 - report, don't abort the batch
                rows.append({"method": meth.label, "pdb_id": pid, "status": "error",
                             "error": str(e)[:150]})
        if progress:
            progress(f"{meth.label}: aligned {n_ok}/{len(folders)}")
    return pd.DataFrame(rows)


def aligned_path(out_root: str, method_label: str, pid: str) -> str | None:
    meth = method_by_label(method_label)
    p = os.path.join(out_root, meth.key, f"{pid}_aligned.pdb")
    return p if os.path.exists(p) else None


# ------------------------------------------------------- pseudosequence map ---
def map_pseudosequences(ds: Dataset, aligned_root: str, progress=None):
    """Map the NetMHCpan-4.1 pseudosequence onto GT and every aligned prediction.

    Returns (gt_df, pred_df, gt_positions, pred_positions). Positions dicts give
    per-structure query residue numbers for the 34 interface positions.
    """
    from histo_pseudosequence.core import PseudosequenceMapper

    mapper = PseudosequenceMapper()
    gt_dir = ds.gt_aligned_dir()

    gt_rows = []
    gt_positions: dict[str, list] = {}
    for pid in sorted(ds.gt_files):
        gf = ds.gt_files[pid][0][1]
        try:
            r = mapper.map(os.path.join(gt_dir, gf))
        except Exception:  # noqa: BLE001
            continue
        gt_positions[pid] = [p.query_resseq for p in r.positions]
        gt_rows.append({"pdb_id": pid, "gt_file": gf, "chain": r.chain,
                        "n_mapped": r.n_mapped, "pseudosequence": r.pseudosequence})
    if progress:
        progress(f"pseudosequence: mapped {len(gt_rows)} GT structures")

    pred_rows = []
    pred_positions: dict[str, list] = {}
    for meth in METHODS:
        n = 0
        for pid in sorted(ds.method_folders.get(meth.key, {})):
            p = aligned_path(aligned_root, meth.label, pid)
            if not p:
                continue
            try:
                r = mapper.map(p)
            except Exception:  # noqa: BLE001
                continue
            pred_positions[f"{meth.label}|{pid}"] = [pp.query_resseq for pp in r.positions]
            gt_pseudo = next((g["pseudosequence"] for g in gt_rows if g["pdb_id"] == pid), None)
            pred_rows.append({"method": meth.label, "pdb_id": pid, "n_mapped": r.n_mapped,
                              "pseudosequence": r.pseudosequence,
                              "matches_gt": (gt_pseudo is not None and r.pseudosequence == gt_pseudo)})
            n += 1
        if progress:
            progress(f"pseudosequence: mapped {n} {meth.label} predictions")
    return pd.DataFrame(gt_rows), pd.DataFrame(pred_rows), gt_positions, pred_positions


# --------------------------------------------------------- peptide scoring ---
def _best_gt_peptide(ds: Dataset, pid: str, pred_pep: dict):
    """GT peptide-atom dict of the copy with lowest whole-peptide Cα RMSD."""
    gt_dir = ds.gt_aligned_dir()
    best = None
    for _copy, gf in ds.gt_files.get(pid, []):
        gm = load_model(os.path.join(gt_dir, gf))
        gp = chain_atoms(gm, PEPTIDE_CHAIN)
        r_ca, _ = rmsd_over(gp, pred_pep, atomset={"CA"})
        if np.isnan(r_ca):
            continue
        if best is None or r_ca < best[0]:
            best = (r_ca, gp, gf)
    return best


def score_peptide_rmsd(ds: Dataset, aligned_root: str, progress=None):
    """Peptide RMSD vs best GT copy: whole (CA + all-atom), per-residue, per-atom."""
    whole, perres, peratom = [], [], []
    for meth in METHODS:
        for pid in sorted(ds.method_folders.get(meth.key, {})):
            p = aligned_path(aligned_root, meth.label, pid)
            if not p:
                continue
            pm = load_model(p)
            pep = chain_atoms(pm, PEPTIDE_CHAIN)
            if not pep:
                continue
            best = _best_gt_peptide(ds, pid, pep)
            if best is None:
                continue
            _rca, gp, gf = best
            shared = sorted(set(gp) & set(pep))
            r_ca, n_ca = rmsd_over(gp, pep, atomset={"CA"})
            r_aa, n_aa = rmsd_over(gp, pep)
            whole.append({"method": meth.label, "pdb_id": pid, "gt_file": gf,
                          "peptide_length": len(shared),
                          "rmsd_ca": round(r_ca, 4), "n_ca": n_ca,
                          "rmsd_allatom": round(r_aa, 4), "n_atoms": n_aa})
            for rn in shared:
                cd = ca_dev(gp[rn], pep[rn])
                aar, naa = rmsd_over({rn: gp[rn]}, {rn: pep[rn]})
                perres.append({"method": meth.label, "pdb_id": pid, "resnum": rn,
                               "ca_dev": round(cd, 4) if cd == cd else None,
                               "allatom_rmsd": round(aar, 4) if aar == aar else None})
                for an in sorted(set(gp[rn]) & set(pep[rn])):
                    peratom.append({"method": meth.label, "pdb_id": pid, "resnum": rn,
                                    "atom": an, "dev": round(float(np.linalg.norm(gp[rn][an] - pep[rn][an])), 4)})
        if progress:
            progress(f"peptide RMSD: {meth.label} done")
    return pd.DataFrame(whole), pd.DataFrame(perres), pd.DataFrame(peratom)


def score_peptide_confidence(ds: Dataset, aligned_root: str, progress=None):
    """Per-residue / per-atom peptide pLDDT (b-factor, ESMFold scaled to 0-100)."""
    whole, perres, peratom = [], [], []
    for meth in METHODS:
        for pid in sorted(ds.method_folders.get(meth.key, {})):
            p = aligned_path(aligned_root, meth.label, pid)
            if not p:
                continue
            m = load_model(p)
            if PEPTIDE_CHAIN not in [c.id for c in m]:
                continue
            res_means = []
            for r in m[PEPTIDE_CHAIN]:
                if r.id[0] != " ":
                    continue
                atoms = [(a.name, a.bfactor * meth.plddt_scale) for a in r if a.element != "H"]
                if not atoms:
                    continue
                rm = float(np.mean([b for _, b in atoms]))
                res_means.append(rm)
                perres.append({"method": meth.label, "pdb_id": pid, "resnum": r.id[1], "plddt": round(rm, 3)})
                for an, b in atoms:
                    peratom.append({"method": meth.label, "pdb_id": pid, "resnum": r.id[1],
                                    "atom": an, "plddt": round(float(b), 3)})
            if res_means:
                whole.append({"method": meth.label, "pdb_id": pid,
                              "peptide_plddt_mean": round(float(np.mean(res_means)), 3),
                              "peptide_plddt_min": round(float(np.min(res_means)), 3),
                              "n_res": len(res_means)})
        if progress:
            progress(f"peptide confidence: {meth.label} done")
    return pd.DataFrame(whole), pd.DataFrame(perres), pd.DataFrame(peratom)


# ------------------------------------------------------- interface scoring ---
def score_interface_rmsd(ds: Dataset, aligned_root: str, gt_positions: dict, progress=None):
    """MHC interface RMSD at pseudosequence positions: mainchain vs sidechain."""
    gt_dir = ds.gt_aligned_dir()
    whole, perpos = [], []
    for meth in METHODS:
        for pid in sorted(ds.method_folders.get(meth.key, {})):
            if pid not in gt_positions:
                continue
            positions = [p for p in gt_positions[pid] if p is not None]
            p = aligned_path(aligned_root, meth.label, pid)
            if not p:
                continue
            prA = chain_atoms(load_model(p), MHC_CHAIN)
            if not prA:
                continue
            best = None
            for _copy, gf in ds.gt_files[pid]:
                gtA = chain_atoms(load_model(os.path.join(gt_dir, gf)), MHC_CHAIN)
                r_mc, _ = rmsd_over(gtA, prA, positions, atomset=MAINCHAIN_ATOMS)
                if np.isnan(r_mc):
                    continue
                if best is None or r_mc < best[0]:
                    best = (r_mc, gtA, gf)
            if best is None:
                continue
            _rmc, gtA, gf = best
            r_mc, n_mc = rmsd_over(gtA, prA, positions, atomset=MAINCHAIN_ATOMS)
            r_sc, n_sc = rmsd_over(gtA, prA, positions, exclude=MAINCHAIN_ATOMS)
            r_ca, _ = rmsd_over(gtA, prA, positions, atomset={"CA"})
            r_all, _ = rmsd_over(gtA, prA, positions)
            whole.append({"method": meth.label, "pdb_id": pid, "gt_file": gf,
                          "n_positions": len(positions),
                          "rmsd_mainchain": round(r_mc, 4), "n_mc_atoms": n_mc,
                          "rmsd_sidechain": round(r_sc, 4) if r_sc == r_sc else None, "n_sc_atoms": n_sc,
                          "rmsd_ca": round(r_ca, 4), "rmsd_allheavy": round(r_all, 4)})
            for rn in positions:
                if rn not in gtA or rn not in prA:
                    continue
                mc, _ = rmsd_over(gtA, prA, [rn], atomset=MAINCHAIN_ATOMS)
                sc, _ = rmsd_over(gtA, prA, [rn], exclude=MAINCHAIN_ATOMS)
                perpos.append({"method": meth.label, "pdb_id": pid, "resnum": rn,
                               "rmsd_mainchain": round(mc, 4) if mc == mc else None,
                               "rmsd_sidechain": round(sc, 4) if sc == sc else None})
        if progress:
            progress(f"interface RMSD: {meth.label} done")
    return pd.DataFrame(whole), pd.DataFrame(perpos)


def score_interface_confidence(ds: Dataset, aligned_root: str, pred_positions: dict, progress=None):
    """MHC interface pLDDT at pseudosequence positions: mainchain vs sidechain."""
    whole, perpos = [], []
    for meth in METHODS:
        for pid in sorted(ds.method_folders.get(meth.key, {})):
            key = f"{meth.label}|{pid}"
            if key not in pred_positions:
                continue
            positions = [p for p in pred_positions[key] if p is not None]
            p = aligned_path(aligned_root, meth.label, pid)
            if not p:
                continue
            m = load_model(p)
            if MHC_CHAIN not in [c.id for c in m]:
                continue
            A = {r.id[1]: r for r in m[MHC_CHAIN] if r.id[0] == " "}
            mc_all, sc_all = [], []
            for rn in positions:
                if rn not in A:
                    continue
                r = A[rn]
                mc = [a.bfactor * meth.plddt_scale for a in r if a.element != "H" and a.name in MAINCHAIN_ATOMS]
                scv = [a.bfactor * meth.plddt_scale for a in r if a.element != "H" and a.name not in MAINCHAIN_ATOMS]
                mc_all += mc
                sc_all += scv
                perpos.append({"method": meth.label, "pdb_id": pid, "resnum": rn,
                               "plddt_mainchain": round(float(np.mean(mc)), 3) if mc else None,
                               "plddt_sidechain": round(float(np.mean(scv)), 3) if scv else None})
            if mc_all or sc_all:
                whole.append({"method": meth.label, "pdb_id": pid, "n_positions": len(positions),
                              "iface_plddt_mainchain": round(float(np.mean(mc_all)), 3) if mc_all else None,
                              "iface_plddt_sidechain": round(float(np.mean(sc_all)), 3) if sc_all else None,
                              "iface_plddt_all": round(float(np.mean(mc_all + sc_all)), 3)})
        if progress:
            progress(f"interface confidence: {meth.label} done")
    return pd.DataFrame(whole), pd.DataFrame(perpos)


# --------------------------------------------------------------- gap + calib ---
def gap_analysis(ds: Dataset):
    """One DataFrame per method of GT complexes it did NOT predict, in the
    pdb_code,locus,allele_slug,peptide_sequence,resolution schema."""
    out: dict[str, pd.DataFrame] = {}
    cols = ["pdb_code", "locus", "allele_slug", "peptide_sequence", "resolution"]
    for meth in METHODS:
        predicted = set(ds.method_folders.get(meth.key, {}))
        missing = sorted(ds.gt_ids - predicted)
        rows = []
        for pid in missing:
            al, pep = ds.pdb_meta.get(pid, (None, None))
            if al is None:
                rows.append({"pdb_code": pid, "locus": "", "allele_slug": "",
                             "peptide_sequence": "", "resolution": 0.0})
            else:
                locus = "_".join(al.split("_")[:2])
                rows.append({"pdb_code": pid, "locus": locus, "allele_slug": al,
                             "peptide_sequence": pep, "resolution": 0.0})
        out[meth.slug] = pd.DataFrame(rows, columns=cols)
    return out


def calibration(conf_whole: pd.DataFrame, rmsd_whole: pd.DataFrame,
                conf_col: str, rmsd_col: str) -> pd.DataFrame:
    """Spearman correlation of a confidence column vs an RMSD column per method
    (negative = well-calibrated: higher confidence -> lower error)."""
    from scipy.stats import spearmanr

    mg = conf_whole.merge(rmsd_whole, on=["method", "pdb_id"])
    rows = []
    for m in METHOD_ORDER:
        s = mg[mg.method == m].dropna(subset=[conf_col, rmsd_col])
        if len(s) < 3:
            rows.append({"method": m, "n": len(s), "spearman": None, "p_value": None})
            continue
        rho, p = spearmanr(s[conf_col], s[rmsd_col])
        rows.append({"method": m, "n": len(s), "spearman": round(float(rho), 3), "p_value": float(p)})
    return pd.DataFrame(rows)
