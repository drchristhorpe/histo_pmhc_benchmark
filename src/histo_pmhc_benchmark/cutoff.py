"""Training-cutoff analysis: split each method's predictions into structures
released before vs after that method's PDB training-data cutoff, and compare
accuracy and confidence across the boundary.

Release dates come from a local metadata CSV when present, with an optional
RCSB/PDBe fallback for PDB codes the CSV doesn't cover.
"""

from __future__ import annotations

import json
import os
import urllib.request

import numpy as np
import pandas as pd

from histo_pmhc_benchmark.config import METHODS, METHOD_ORDER, method_by_label

# Default location of the release-date CSV inside a dataset dir (pdb_code,release_date).
_RELEASE_CSV = os.path.join("metadata", "alphafold_date_cutoff_pdb_codes.csv")


def _fetch_release_date(pdb_id: str, timeout: float = 15.0) -> pd.Timestamp | None:
    """Fetch a PDB entry's initial release date from RCSB, then PDBe. None if
    unavailable (e.g. obsoleted entry)."""
    pid = pdb_id.lower()
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pid.upper()}"
        with urllib.request.urlopen(url, timeout=timeout) as r:
            d = json.load(r)
        rel = d.get("rcsb_accession_info", {}).get("initial_release_date")
        if rel:
            return pd.Timestamp(rel).tz_localize(None)
    except Exception:  # noqa: BLE001
        pass
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pid}"
        with urllib.request.urlopen(url, timeout=timeout) as r:
            d = json.load(r)
        rel = d[pid][0].get("release_date")
        if rel:
            return pd.Timestamp(rel)
    except Exception:  # noqa: BLE001
        pass
    return None


def load_release_dates(base: str, pdb_ids=None, fetch_missing: bool = True,
                       progress=None) -> dict[str, pd.Timestamp]:
    """Release-date map for pdb_ids: local metadata CSV first, RCSB/PDBe for the rest.

    fetch_missing=False stays fully offline (missing ids are simply absent, and
    downstream classification marks them 'Unknown').
    """
    reldate: dict[str, pd.Timestamp] = {}
    csv_path = os.path.join(base, _RELEASE_CSV)
    if os.path.exists(csv_path):
        rd = pd.read_csv(csv_path)
        rd["pdb_code"] = rd["pdb_code"].str.lower()
        for pid, d in zip(rd["pdb_code"], pd.to_datetime(rd["release_date"])):
            reldate[pid] = d
    if pdb_ids is not None and fetch_missing:
        missing = sorted({p.lower() for p in pdb_ids} - set(reldate))
        if missing and progress:
            progress(f"cutoff: fetching {len(missing)} missing release dates from RCSB/PDBe")
        for pid in missing:
            got = _fetch_release_date(pid)
            if got is not None:
                reldate[pid] = got
    return reldate


def classify(df: pd.DataFrame, reldate: dict[str, pd.Timestamp]) -> pd.DataFrame:
    """Add release_date, cutoff and period (Before/After/Unknown) columns keyed on
    each row's method + pdb_id."""
    cutoffs = {m.label: pd.Timestamp(m.training_cutoff) for m in METHODS}
    out = df.copy()
    out["release_date"] = out["pdb_id"].map(lambda p: reldate.get(str(p).lower()))
    out["cutoff"] = out["method"].map(cutoffs)
    out["period"] = np.where(
        out["release_date"].notna(),
        np.where(out["release_date"] < out["cutoff"], "Before", "After"),
        "Unknown",
    )
    return out


def summarize(df: pd.DataFrame, value_col: str, metric_label: str,
              min_n: int = 5) -> pd.DataFrame:
    """Before-vs-after summary for one metric: per-method medians, ratio, and a
    Mann-Whitney U two-sided p-value (computed only when both arms have >=min_n)."""
    from scipy.stats import mannwhitneyu

    rows = []
    for m in METHOD_ORDER:
        sub = df[(df["method"] == m) & (df["period"].isin(["Before", "After"]))].dropna(subset=[value_col])
        b = sub[sub["period"] == "Before"][value_col]
        a = sub[sub["period"] == "After"][value_col]
        rec = {"metric": metric_label, "method": m, "n_before": len(b), "n_after": len(a),
               "median_before": round(float(b.median()), 3) if len(b) else None,
               "median_after": round(float(a.median()), 3) if len(a) else None,
               "ratio_after_over_before": None, "mwu_p": None}
        if len(b) >= min_n and len(a) >= min_n:
            if b.median():
                rec["ratio_after_over_before"] = round(float(a.median() / b.median()), 3)
            _, p = mannwhitneyu(b, a, alternative="two-sided")
            rec["mwu_p"] = round(float(p), 4)
        rows.append(rec)
    return pd.DataFrame(rows)


def run_cutoff_analysis(base: str, tables: dict[str, pd.DataFrame],
                        fetch_missing: bool = True, progress=None) -> dict[str, pd.DataFrame]:
    """Full cutoff analysis over already-scored tables.

    tables: dict with keys peptide_rmsd_whole, interface_rmsd_whole,
            peptide_confidence_whole, interface_confidence_whole, and optionally the
            per-residue/per-position tables (peptide_confidence_per_residue,
            interface_confidence_per_position, peptide_rmsd_per_residue,
            interface_rmsd_per_position) — classified copies of any present are returned.
    Returns {classification, accuracy_summary, confidence_summary, *_classified, *_cutoff}.
    """
    log = progress or (lambda *_: None)
    pep = tables["peptide_rmsd_whole"]
    reldate = load_release_dates(base, pdb_ids=set(pep["pdb_id"]), fetch_missing=fetch_missing, progress=log)
    log(f"cutoff: {len(reldate)} release dates resolved")

    pepC = classify(pep, reldate)
    ifaceC = classify(tables["interface_rmsd_whole"], reldate)
    pconfC = classify(tables["peptide_confidence_whole"], reldate)
    iconfC = classify(tables["interface_confidence_whole"], reldate)

    accuracy = pd.concat([
        summarize(pepC, "rmsd_ca", "Peptide Cα RMSD (Å)"),
        summarize(ifaceC, "rmsd_mainchain", "Interface mainchain RMSD (Å)"),
        summarize(ifaceC, "rmsd_sidechain", "Interface sidechain RMSD (Å)"),
    ], ignore_index=True)
    confidence = pd.concat([
        summarize(pconfC, "peptide_plddt_mean", "Peptide pLDDT"),
        summarize(iconfC, "iface_plddt_mainchain", "Interface mainchain pLDDT"),
        summarize(iconfC, "iface_plddt_sidechain", "Interface sidechain pLDDT"),
    ], ignore_index=True)

    out = {
        "classification": pepC[["method", "pdb_id", "release_date", "cutoff", "period"]],
        "accuracy_summary": accuracy,
        "confidence_summary": confidence,
        # classified whole-confidence tables, so figures reuse them without re-fetching
        "peptide_confidence_classified": pconfC,
        "interface_confidence_classified": iconfC,
    }
    # classify any per-residue / per-position tables that were supplied
    for key in ("peptide_confidence_per_residue", "interface_confidence_per_position",
                "peptide_rmsd_per_residue", "interface_rmsd_per_position"):
        if key in tables:
            out[key + "_cutoff"] = classify(tables[key], reldate)
    return out


def add_pseudo_position_index(df: pd.DataFrame) -> pd.DataFrame:
    """Add a canonical pseudosequence-position rank (pos_idx, 1..34) so per-structure
    residue numbers map onto a common x-axis for per-position figures."""
    out = df.copy()
    out["pos_idx"] = out.groupby(["method", "pdb_id"])["resnum"].rank(method="dense").astype(int)
    return out
