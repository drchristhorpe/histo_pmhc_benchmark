"""End-to-end orchestration: run all stages against a dataset and write outputs."""

from __future__ import annotations

import json
import os

import pandas as pd

from histo_pmhc_benchmark.dataset import build_dataset
from histo_pmhc_benchmark import stages
from histo_pmhc_benchmark import cutoff as cutoff_mod
from histo_pmhc_benchmark.figures import (
    summary_figure, confidence_cutoff_boxplots, per_residue_cutoff_figure,
)
from histo_pmhc_benchmark.config import method_by_label


def _write(df: pd.DataFrame, out_dir: str, name: str) -> str:
    path = os.path.join(out_dir, name)
    df.to_csv(path, index=False)
    return path


def _per_residue_cutoff_figures(cutoff_out, peptide_whole, out_dir, written, log):
    """Render one per-residue before/after box-plot figure per method that has an
    after-cutoff set. Peptide panels use the dominant peptide length; interface
    panels use the canonical pseudosequence-position rank."""
    from histo_pmhc_benchmark import cutoff as _cut

    pc = cutoff_out.get("peptide_confidence_per_residue_cutoff")
    pr = cutoff_out.get("peptide_rmsd_per_residue_cutoff")
    ic = cutoff_out.get("interface_confidence_per_position_cutoff")
    ir = cutoff_out.get("interface_rmsd_per_position_cutoff")
    if any(x is None for x in (pc, pr, ic, ir)):
        log("per-residue cutoff figures skipped: missing per-residue tables")
        return
    ic = _cut.add_pseudo_position_index(ic)
    ir = _cut.add_pseudo_position_index(ir)

    # dominant peptide length among scored complexes (for a clean common x-axis)
    plen = int(peptide_whole["peptide_length"].mode().iloc[0]) if "peptide_length" in peptide_whole else 9
    pep9 = set(peptide_whole[peptide_whole["peptide_length"] == plen]["pdb_id"]) \
        if "peptide_length" in peptide_whole else set(peptide_whole["pdb_id"])

    classification = cutoff_out["classification"]
    for meth in sorted(classification["method"].unique()):
        n_after = int(((classification["method"] == meth) & (classification["period"] == "After")).sum())
        n_before = int(((classification["method"] == meth) & (classification["period"] == "Before")).sum())
        if n_after < 5:
            continue  # nothing meaningful to compare
        slug = method_by_label(meth).slug
        out_path = os.path.join(out_dir, f"pmhc_perres_cutoff_{slug}.png")
        per_residue_cutoff_figure(
            meth,
            pc[(pc.method == meth) & (pc.pdb_id.isin(pep9)) & (pc.resnum <= plen)],
            pr[(pr.method == meth) & (pr.pdb_id.isin(pep9)) & (pr.resnum <= plen)],
            ic[ic.method == meth], ir[ir.method == meth], out_path,
            peptide_length=plen, n_before=n_before, n_after=n_after)
        written[f"perres_cutoff_{slug}"] = out_path
        log(f"per-residue cutoff figure: {meth}")


def run_pipeline(base: str, out_dir: str, aligned_root: str | None = None,
                 do_align: bool = True, do_figure: bool = True,
                 do_cutoff: bool = True, fetch_missing_dates: bool = True,
                 progress=None) -> dict:
    """Run the full benchmark. Returns a manifest of written files + headline stats.

    base:        dataset dir (predictions/, ground_truth/aligned/, metadata/).
    out_dir:     where CSVs and the figure are written.
    aligned_root: where aligned PDBs live/are written (default <out_dir>/aligned).
    do_align:    if False, reuse existing aligned PDBs under aligned_root.
    do_cutoff:   run the before/after training-cutoff analysis (needs release dates).
    fetch_missing_dates: allow RCSB/PDBe lookup for PDB codes the local metadata
                 CSV doesn't cover (set False to stay fully offline).
    """
    os.makedirs(out_dir, exist_ok=True)
    aligned_root = aligned_root or os.path.join(out_dir, "aligned")
    log = progress or (lambda *_: None)

    ds = build_dataset(base)
    log(f"dataset: {len(ds.gt_ids)} GT complexes, "
        f"{sum(len(v) for v in ds.method_folders.values())} prediction folders")
    written: dict[str, str] = {}

    # 1. align
    if do_align:
        rep = stages.align_predictions(ds, aligned_root, progress=log)
        written["alignment_report"] = _write(rep, out_dir, "pmhc_alignment_report.csv")

    # 2. pseudosequence
    gt_ps, pred_ps, gt_pos, pred_pos = stages.map_pseudosequences(ds, aligned_root, progress=log)
    written["gt_pseudosequence"] = _write(gt_ps, out_dir, "pmhc_gt_pseudosequence.csv")
    written["prediction_pseudosequence"] = _write(pred_ps, out_dir, "pmhc_prediction_pseudosequence.csv")

    # 3. peptide RMSD + confidence
    pw, pr, pa = stages.score_peptide_rmsd(ds, aligned_root, progress=log)
    written["peptide_rmsd_whole"] = _write(pw, out_dir, "pmhc_peptide_rmsd_whole.csv")
    written["peptide_rmsd_per_residue"] = _write(pr, out_dir, "pmhc_peptide_rmsd_per_residue.csv")
    written["peptide_rmsd_per_atom"] = _write(pa, out_dir, "pmhc_peptide_rmsd_per_atom.csv")
    cw, cr, ca = stages.score_peptide_confidence(ds, aligned_root, progress=log)
    written["peptide_confidence_whole"] = _write(cw, out_dir, "pmhc_peptide_confidence_whole.csv")
    written["peptide_confidence_per_residue"] = _write(cr, out_dir, "pmhc_peptide_confidence_per_residue.csv")
    written["peptide_confidence_per_atom"] = _write(ca, out_dir, "pmhc_peptide_confidence_per_atom.csv")
    pep_cal = stages.calibration(cw, pw, "peptide_plddt_mean", "rmsd_ca")
    written["peptide_confidence_calibration"] = _write(pep_cal, out_dir, "pmhc_peptide_confidence_calibration.csv")

    # 4. interface RMSD + confidence
    iw, ip = stages.score_interface_rmsd(ds, aligned_root, gt_pos, progress=log)
    written["interface_rmsd_whole"] = _write(iw, out_dir, "pmhc_interface_rmsd_whole.csv")
    written["interface_rmsd_per_position"] = _write(ip, out_dir, "pmhc_interface_rmsd_per_position.csv")
    icw, icp = stages.score_interface_confidence(ds, aligned_root, pred_pos, progress=log)
    written["interface_confidence_whole"] = _write(icw, out_dir, "pmhc_interface_confidence_whole.csv")
    written["interface_confidence_per_position"] = _write(icp, out_dir, "pmhc_interface_confidence_per_position.csv")
    ifc_cal = stages.calibration(icw, iw, "iface_plddt_sidechain", "rmsd_sidechain")
    written["interface_confidence_calibration"] = _write(ifc_cal, out_dir, "pmhc_interface_confidence_calibration.csv")

    # 5. gap analysis
    for slug, df in stages.gap_analysis(ds).items():
        written[f"missing_{slug}"] = _write(df, out_dir, f"missing_{slug}.csv")

    # 6. training-cutoff analysis (before vs after each method's PDB cutoff)
    cutoff_out = None
    if do_cutoff:
        try:
            cutoff_out = cutoff_mod.run_cutoff_analysis(
                base,
                {"peptide_rmsd_whole": pw, "interface_rmsd_whole": iw,
                 "peptide_confidence_whole": cw, "interface_confidence_whole": icw,
                 "peptide_confidence_per_residue": cr, "interface_confidence_per_position": icp,
                 "peptide_rmsd_per_residue": pr, "interface_rmsd_per_position": ip},
                fetch_missing=fetch_missing_dates, progress=log)
            written["cutoff_classification"] = _write(cutoff_out["classification"], out_dir, "pmhc_cutoff_classification.csv")
            written["cutoff_accuracy_summary"] = _write(cutoff_out["accuracy_summary"], out_dir, "pmhc_cutoff_accuracy_summary.csv")
            written["cutoff_confidence_summary"] = _write(cutoff_out["confidence_summary"], out_dir, "pmhc_cutoff_confidence_summary.csv")
            if "peptide_confidence_per_residue_cutoff" in cutoff_out:
                written["peptide_confidence_per_residue_cutoff"] = _write(
                    cutoff_out["peptide_confidence_per_residue_cutoff"], out_dir,
                    "pmhc_peptide_confidence_per_residue_cutoff.csv")
            if "interface_confidence_per_position_cutoff" in cutoff_out:
                written["interface_confidence_per_position_cutoff"] = _write(
                    cutoff_out["interface_confidence_per_position_cutoff"], out_dir,
                    "pmhc_interface_confidence_per_position_cutoff.csv")
        except Exception as e:  # noqa: BLE001
            log(f"cutoff analysis skipped: {e}")

    # 7. figures
    if do_figure and len(pw):
        try:
            summary_figure(pw, pr, iw, pep_cal, ifc_cal,
                           os.path.join(out_dir, "pmhc_summary_figure.png"))
            written["summary_figure"] = os.path.join(out_dir, "pmhc_summary_figure.png")
        except Exception as e:  # noqa: BLE001
            log(f"summary figure skipped: {e}")
        if cutoff_out is not None:
            try:
                confidence_cutoff_boxplots(
                    cutoff_out["peptide_confidence_classified"],
                    cutoff_out["interface_confidence_classified"],
                    os.path.join(out_dir, "pmhc_confidence_cutoff_boxplots.png"))
                written["confidence_cutoff_boxplots"] = os.path.join(out_dir, "pmhc_confidence_cutoff_boxplots.png")
            except Exception as e:  # noqa: BLE001
                log(f"cutoff box plots skipped: {e}")
            # per-residue before/after box plots, one figure per method with an after set
            try:
                _per_residue_cutoff_figures(cutoff_out, pw, out_dir, written, log)
            except Exception as e:  # noqa: BLE001
                log(f"per-residue cutoff figures skipped: {e}")

    # headline stats
    headline = {}
    for m in pw.method.unique():
        s = pw[pw.method == m]; i = iw[iw.method == m]
        headline[m] = {
            "n": int(len(s)),
            "peptide_ca_median": round(float(s.rmsd_ca.median()), 3),
            "interface_mainchain_median": round(float(i.rmsd_mainchain.median()), 3) if len(i) else None,
            "interface_sidechain_median": round(float(i.rmsd_sidechain.median()), 3) if len(i) else None,
        }
    manifest = {"out_dir": out_dir, "aligned_root": aligned_root,
                "n_gt": len(ds.gt_ids), "files": written, "headline": headline}
    json.dump(manifest, open(os.path.join(out_dir, "pmhc_run_manifest.json"), "w"), indent=2)
    written["manifest"] = os.path.join(out_dir, "pmhc_run_manifest.json")
    return manifest
