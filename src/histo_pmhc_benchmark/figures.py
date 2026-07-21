"""Summary figure for the pMHC benchmark (4 panels)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from histo_pmhc_benchmark.config import METHOD_COLORS, METHOD_ORDER


def summary_figure(peptide_whole: pd.DataFrame, peptide_perres: pd.DataFrame,
                   interface_whole: pd.DataFrame, pep_calib: pd.DataFrame,
                   iface_calib: pd.DataFrame, out_path: str):
    """Render the 4-panel benchmark summary and save to out_path (PNG)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order = [m for m in METHOD_ORDER if m in set(peptide_whole.method)]
    mc = {m: METHOD_COLORS[m] for m in order}
    grey = "#888888"

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axA, axB, axC, axD = axes.ravel()

    # A: peptide CA RMSD violins
    data = [peptide_whole[peptide_whole.method == m].rmsd_ca.dropna().values for m in order]
    parts = axA.violinplot(data, showextrema=False, widths=0.8)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(mc[order[i]]); pc.set_alpha(0.55)
        pc.set_edgecolor(mc[order[i]]); pc.set_linewidth(1.0)
    for i, d in enumerate(data):
        if len(d):
            md = float(np.median(d))
            axA.hlines(md, i + 0.72, i + 1.28, color="black", lw=2, zorder=5)
            axA.text(i + 1, md + 0.06, f"{md:.2f}", ha="center", va="bottom", fontsize=6, fontweight="bold")
    axA.axhline(1.0, color=grey, ls=":", lw=0.8)
    axA.set_xticks(range(1, len(order) + 1))
    axA.set_xticklabels([o.replace("-", "-\n") if len(o) > 9 else o for o in order], fontsize=6)
    axA.set_ylabel("Peptide Cα RMSD (Å)"); axA.set_ylim(0, 3.0)
    axA.set_title("Peptide backbone accuracy (whole peptide)", fontsize=8, loc="left")

    # B: per-residue CA deviation profile for 9-mers
    nmer9 = set(peptide_whole[peptide_whole.peptide_length == 9].pdb_id)
    p9 = peptide_perres[peptide_perres.pdb_id.isin(nmer9) & (peptide_perres.resnum <= 9)]
    for m in order:
        prof = p9[p9.method == m].groupby("resnum")["ca_dev"].median()
        if len(prof):
            axB.plot(prof.index, prof.values, "-o", color=mc[m], lw=1.5, ms=4, label=m)
    axB.set_xlabel("Peptide residue position (9-mers)")
    axB.set_ylabel("Median Cα deviation (Å)")
    axB.set_xticks(range(1, 10))
    axB.set_title("Central bulge is hardest; anchors are tight", fontsize=8, loc="left")
    axB.legend(frameon=False, fontsize=6, loc="upper right")

    # C: interface mainchain vs sidechain
    x = np.arange(len(order)); w = 0.36
    mcv = [interface_whole[interface_whole.method == m].rmsd_mainchain.median() for m in order]
    scv = [interface_whole[interface_whole.method == m].rmsd_sidechain.median() for m in order]
    axC.bar(x - w / 2, mcv, w, color=[mc[m] for m in order], alpha=0.95)
    axC.bar(x + w / 2, scv, w, color=[mc[m] for m in order], alpha=0.45, hatch="///")
    axC.set_xticks(x)
    axC.set_xticklabels([o.replace("-", "-\n") if len(o) > 9 else o for o in order], fontsize=6)
    axC.set_ylabel("Interface RMSD (Å)")
    axC.set_title("MHC binding groove: backbone vs sidechain", fontsize=8, loc="left")
    axC.legend(handles=[Patch(facecolor=grey, alpha=0.95, label="Mainchain (N,Cα,C,O)"),
                        Patch(facecolor=grey, alpha=0.45, hatch="///", label="Sidechain")],
               frameon=False, fontsize=6, loc="upper left")

    # D: calibration magnitude (|rho|) for peptide + interface sidechain
    def rho(df, m):
        v = df[df.method == m]["spearman"]
        return -float(v.values[0]) if len(v) and v.values[0] is not None else 0.0
    pep_rho = [rho(pep_calib, m) for m in order]
    ifc_rho = [rho(iface_calib, m) for m in order]
    axD.bar(x - w / 2, pep_rho, w, color=[mc[m] for m in order], alpha=0.95)
    axD.bar(x + w / 2, ifc_rho, w, color=[mc[m] for m in order], alpha=0.45, hatch="///")
    axD.set_xticks(x)
    axD.set_xticklabels([o.replace("-", "-\n") if len(o) > 9 else o for o in order], fontsize=6)
    axD.set_ylabel("Calibration |Spearman ρ|\n(pLDDT vs RMSD)")
    axD.set_title("Confidence↔accuracy calibration", fontsize=8, loc="left")
    axD.legend(handles=[Patch(facecolor=grey, alpha=0.95, label="Peptide (Cα)"),
                        Patch(facecolor=grey, alpha=0.45, hatch="///", label="Interface sidechain")],
               frameon=False, fontsize=6, loc="upper right")
    axD.set_ylim(0, max([*pep_rho, *ifc_rho, 0.1]) * 1.25)

    for ax, l in zip([axA, axB, axC, axD], "abcd"):
        ax.text(-0.08, 1.04, l, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    n = len(peptide_whole)
    fig.suptitle(f"pMHC class I structure prediction: peptide and MHC-groove accuracy "
                 f"(n={n} predictions, {len(order)} methods)", fontsize=9, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    return out_path
