"""End-to-end orchestration: run all stages against a dataset and write outputs."""

from __future__ import annotations

import json
import os

import pandas as pd

from histo_pmhc_benchmark.dataset import build_dataset
from histo_pmhc_benchmark import stages
from histo_pmhc_benchmark.figures import summary_figure


def _write(df: pd.DataFrame, out_dir: str, name: str) -> str:
    path = os.path.join(out_dir, name)
    df.to_csv(path, index=False)
    return path


def run_pipeline(base: str, out_dir: str, aligned_root: str | None = None,
                 do_align: bool = True, do_figure: bool = True, progress=None) -> dict:
    """Run the full benchmark. Returns a manifest of written files + headline stats.

    base:        dataset dir (predictions/, ground_truth/aligned/, metadata/).
    out_dir:     where CSVs and the figure are written.
    aligned_root: where aligned PDBs live/are written (default <out_dir>/aligned).
    do_align:    if False, reuse existing aligned PDBs under aligned_root.
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

    # 6. figure
    if do_figure and len(pw):
        try:
            summary_figure(pw, pr, iw, pep_cal, ifc_cal,
                           os.path.join(out_dir, "pmhc_summary_figure.png"))
            written["summary_figure"] = os.path.join(out_dir, "pmhc_summary_figure.png")
        except Exception as e:  # noqa: BLE001
            log(f"figure skipped: {e}")

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
