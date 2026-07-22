"""Summary figures for the pMHC benchmark."""

from __future__ import annotations

import numpy as np
import pandas as pd

from histo_pmhc_benchmark.config import METHOD_COLORS, METHOD_ORDER

# before = neutral grey, after = blue — the project's training-cutoff convention.
PERIOD_COLORS = {"Before": "#666666", "After": "#648fff"}


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


def _box_pair(ax, df, valcol, positions, poscol, ylabel, title, min_n=5, xticklabels=None):
    """Paired before/after box plots at each position on one axis. `positions` are the
    values matched in `poscol`; `xticklabels` (optional) relabels the ticks — e.g. the
    canonical MHC residue numbers instead of a 1..34 rank."""
    for i, pos in enumerate(positions):
        for per, dx in [("Before", -0.19), ("After", 0.19)]:
            vals = df[(df[poscol] == pos) & (df["period"] == per)][valcol].dropna().values
            if len(vals) >= min_n:
                bp = ax.boxplot([vals], positions=[i + dx], widths=0.30, patch_artist=True,
                                showfliers=False, medianprops=dict(color="black", lw=1.1),
                                whiskerprops=dict(color=PERIOD_COLORS[per], lw=0.7),
                                capprops=dict(color=PERIOD_COLORS[per], lw=0.7))
                for b in bp["boxes"]:
                    b.set_facecolor(PERIOD_COLORS[per]); b.set_alpha(0.55)
                    b.set_edgecolor(PERIOD_COLORS[per]); b.set_linewidth(0.7)
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(xticklabels if xticklabels is not None else positions, fontsize=6)
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontsize=8)
    ax.margins(x=0.02)


def per_residue_cutoff_figure(method: str, peptide_conf: pd.DataFrame, peptide_rmsd: pd.DataFrame,
                              iface_conf: pd.DataFrame, iface_rmsd: pd.DataFrame, out_path: str,
                              peptide_length: int = 9, iface_residue_numbers=None,
                              n_before: int | None = None, n_after: int | None = None):
    """Per-residue before/after box plots for one method: peptide pLDDT + Cα deviation
    (positions 1..peptide_length) and binding-site sidechain pLDDT + RMSD at the
    NetMHCpan-4.1 positions. Each input table must carry a `period` column; the interface
    tables must carry a `pos_idx` (canonical pseudosequence rank). The interface x-axis is
    labelled with the canonical MHC residue numbers (iface_residue_numbers, one per rank).
    Boxes show the variance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    from histo_pseudosequence.reference_data import NETMHCPAN_41_POSITIONS
    if iface_residue_numbers is None:
        iface_residue_numbers = list(NETMHCPAN_41_POSITIONS)
    n_positions = len(iface_residue_numbers)

    pep_pos = list(range(1, peptide_length + 1))
    ifc_pos = list(range(1, n_positions + 1))
    fig, axes = plt.subplots(4, 1, figsize=(13, 13))
    _box_pair(axes[0], peptide_conf, "plddt", pep_pos, "resnum",
              "Peptide pLDDT", "Peptide — per-residue confidence (pLDDT)")
    _box_pair(axes[1], peptide_rmsd, "ca_dev", pep_pos, "resnum",
              "Peptide Cα dev (Å)", "Peptide — per-residue accuracy (Cα deviation)")
    _box_pair(axes[2], iface_conf, "plddt_sidechain", ifc_pos, "pos_idx",
              "Interface sidechain pLDDT", "MHC binding site — per-residue confidence (sidechain pLDDT)",
              xticklabels=iface_residue_numbers)
    _box_pair(axes[3], iface_rmsd, "rmsd_sidechain", ifc_pos, "pos_idx",
              "Interface sidechain RMSD (Å)", "MHC binding site — per-residue accuracy (sidechain RMSD)",
              xticklabels=iface_residue_numbers)
    axes[1].set_xlabel(f"Peptide residue position ({peptide_length}-mers)")
    axes[3].set_xlabel("MHC heavy-chain residue number (NetMHCpan-4.1 contact positions)")
    axes[0].legend(handles=[Patch(facecolor=PERIOD_COLORS["Before"], alpha=0.55, label="Before cutoff"),
                            Patch(facecolor=PERIOD_COLORS["After"], alpha=0.55, label="After cutoff")],
                   frameon=False, fontsize=6, loc="lower left", ncol=2)
    for ax, l in zip(axes, "abcd"):
        ax.text(-0.05, 1.02, l, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    nsuff = ""
    if n_before is not None and n_after is not None:
        nsuff = f" (before n={n_before}, after n={n_after})"
    fig.suptitle(f"{method}: per-residue confidence and accuracy, before vs after training cutoff{nsuff}",
                 fontsize=9, fontweight="bold", y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as _plt
    _plt.close(fig)
    return out_path


def confidence_cutoff_boxplots(pconf_cutoff: pd.DataFrame, iconf_cutoff: pd.DataFrame,
                               out_path: str, min_n: int = 5):
    """Box plots of pLDDT before vs after training cutoff, for three compartments
    (peptide / interface mainchain / interface sidechain), paired per method.

    Box plots make the variance visible: box = IQR, line = median, whiskers = 1.5×IQR,
    fliers = outliers. `pconf_cutoff` / `iconf_cutoff` are the whole-structure
    confidence tables with a `period` column (from cutoff.classify)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order = [m for m in METHOD_ORDER if m in set(pconf_cutoff["method"])]
    compartments = [
        (pconf_cutoff, "peptide_plddt_mean", "Peptide", (74, 101)),
        (iconf_cutoff, "iface_plddt_mainchain", "MHC interface — mainchain", (92, 100.5)),
        (iconf_cutoff, "iface_plddt_sidechain", "MHC interface — sidechain", (92, 100.5)),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(9, 10.5))
    off = 0.19
    for ax, (df, col, title, yl) in zip(axes, compartments):
        for i, m in enumerate(order):
            for per, dx in [("Before", -off), ("After", +off)]:
                vals = df[(df["method"] == m) & (df["period"] == per)][col].dropna().values
                if len(vals) >= min_n:
                    bp = ax.boxplot([vals], positions=[i + dx], widths=0.30, patch_artist=True,
                                    showfliers=True, flierprops=dict(marker="o", ms=2,
                                    markerfacecolor=PERIOD_COLORS[per], markeredgecolor="none", alpha=0.4),
                                    medianprops=dict(color="black", lw=1.4),
                                    whiskerprops=dict(color=PERIOD_COLORS[per], lw=0.9),
                                    capprops=dict(color=PERIOD_COLORS[per], lw=0.9))
                    for box in bp["boxes"]:
                        box.set_facecolor(PERIOD_COLORS[per]); box.set_alpha(0.55)
                        box.set_edgecolor(PERIOD_COLORS[per]); box.set_linewidth(0.9)
                elif len(vals) >= 1:
                    ax.scatter(np.full(len(vals), i + dx), vals, s=8,
                               color=PERIOD_COLORS[per], alpha=0.6, zorder=4)
        ax.set_ylabel("pLDDT"); ax.set_title(title, loc="left")
        ax.set_ylim(*yl); ax.margins(x=0.04)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order if ax is axes[-1] else [""] * len(order))
        # note methods with no after-cutoff data
        for i, m in enumerate(order):
            n_after = len(df[(df["method"] == m) & (df["period"] == "After")][col].dropna())
            if n_after == 0:
                ax.text(i + off, yl[0] + (yl[1] - yl[0]) * 0.04, "no post-\ncutoff data",
                        ha="center", va="bottom", fontsize=5.5, color="#666666", style="italic")
    axes[0].legend(handles=[Patch(facecolor=PERIOD_COLORS["Before"], alpha=0.55, label="Before cutoff"),
                            Patch(facecolor=PERIOD_COLORS["After"], alpha=0.55, label="After cutoff")],
                   frameon=False, fontsize=6, loc="lower right", ncol=2)
    for ax, l in zip(axes, "abc"):
        ax.text(-0.07, 1.03, l, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    fig.suptitle("pMHC prediction confidence (pLDDT) before vs after each method's training cutoff",
                 fontsize=9, fontweight="bold", y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    return out_path
