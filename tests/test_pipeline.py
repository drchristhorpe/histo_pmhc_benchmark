"""Unit tests for the benchmark pipeline components.

These exercise the pure logic (geometry, config, gap schema, calibration) that
does not require the full dataset, using small hand-built inputs and one real
committed fixture pair.
"""

import numpy as np
import pandas as pd

from histo_pmhc_benchmark.config import (
    MAINCHAIN_ATOMS,
    METHODS,
    METHOD_COLORS,
    METHOD_ORDER,
    method_by_key,
    method_by_label,
)
from histo_pmhc_benchmark.geometry import ca_dev, rmsd_over
from histo_pmhc_benchmark.stages import calibration, gap_analysis
from histo_pmhc_benchmark.dataset import Dataset


def test_five_methods_distinct_colors():
    assert len(METHODS) == 5
    assert len(set(METHOD_COLORS.values())) == 5  # no shared hex
    # method identity must not reuse the reserved IBM bin/cutoff hexes
    reserved = {"#785EF0", "#DC267F", "#FE6100", "#FFB000", "#648FFF"}
    assert not (set(METHOD_COLORS.values()) & reserved)


def test_method_lookup():
    assert method_by_key("boltz2").label == "Boltz-2"
    assert method_by_label("HistoFold").key == "histofold"
    assert method_by_key("nope") is None


def test_rmsd_identical_is_zero():
    gt = {1: {"CA": np.array([0.0, 0.0, 0.0]), "CB": np.array([1.0, 0.0, 0.0])}}
    r, n = rmsd_over(gt, gt)
    assert r == 0.0 and n == 2


def test_rmsd_known_offset():
    gt = {1: {"CA": np.array([0.0, 0.0, 0.0])}}
    pred = {1: {"CA": np.array([3.0, 4.0, 0.0])}}
    r, n = rmsd_over(gt, pred, atomset={"CA"})
    assert abs(r - 5.0) < 1e-9 and n == 1


def test_rmsd_mainchain_vs_sidechain_split():
    gt = {1: {"CA": np.zeros(3), "N": np.zeros(3), "CB": np.zeros(3)}}
    pred = {1: {"CA": np.zeros(3), "N": np.zeros(3), "CB": np.array([2.0, 0.0, 0.0])}}
    mc, _ = rmsd_over(gt, pred, atomset=MAINCHAIN_ATOMS)
    sc, _ = rmsd_over(gt, pred, exclude=MAINCHAIN_ATOMS)
    assert mc == 0.0          # backbone perfect
    assert abs(sc - 2.0) < 1e-9  # sidechain offset


def test_ca_dev():
    assert ca_dev({"CA": np.zeros(3)}, {"CA": np.array([0.0, 0.0, 2.0])}) == 2.0
    assert np.isnan(ca_dev({"CB": np.zeros(3)}, {"CA": np.zeros(3)}))


def test_gap_schema():
    ds = Dataset(base="/x",
                 gt_files={"1abc": [(1, "1abc_1.pdb")], "2def": [(1, "2def_1.pdb")]},
                 method_folders={"boltz2": {"1abc": "1abc_x"}},
                 pdb_meta={"2def": ("hla_a_02_01", "SIINFEKL")})
    gaps = gap_analysis(ds)
    bo = gaps["boltz2"]
    assert list(bo.columns) == ["pdb_code", "locus", "allele_slug", "peptide_sequence", "resolution"]
    row = bo[bo.pdb_code == "2def"].iloc[0]
    assert row.locus == "hla_a" and row.allele_slug == "hla_a_02_01"
    assert row.peptide_sequence == "SIINFEKL" and row.resolution == 0.0


def test_calibration_sign():
    # perfectly anti-correlated confidence/rmsd -> spearman -1
    conf = pd.DataFrame({"method": ["Boltz-2"] * 4, "pdb_id": list("abcd"),
                         "plddt": [90, 80, 70, 60]})
    rmsd = pd.DataFrame({"method": ["Boltz-2"] * 4, "pdb_id": list("abcd"),
                         "rmsd": [0.1, 0.2, 0.3, 0.4]})
    cal = calibration(conf, rmsd, "plddt", "rmsd")
    row = cal[cal.method == "Boltz-2"].iloc[0]
    assert abs(row["spearman"] + 1.0) < 1e-9
    # methods with no data are reported with n=0, not dropped
    assert (cal["n"] == 0).sum() == 4
