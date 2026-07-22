"""Tests for the training-cutoff analysis."""

import pandas as pd

from histo_pmhc_benchmark.config import method_by_label
from histo_pmhc_benchmark.cutoff import classify, summarize


def test_histofold_uses_alphafold2_cutoff():
    # HistoFold builds on AlphaFold2 (2018-04-30), NOT AlphaFold3 (2021-09-30)
    assert method_by_label("HistoFold").training_cutoff == "2018-04-30"
    assert method_by_label("AlphaFold3").training_cutoff == "2021-09-30"
    assert method_by_label("Boltz-2").training_cutoff == "2023-06-01"


def test_classify_periods():
    df = pd.DataFrame({"method": ["AlphaFold3", "AlphaFold3", "HistoFold"],
                       "pdb_id": ["1abc", "2def", "3ghi"], "rmsd_ca": [0.4, 0.5, 0.9]})
    reldate = {"1abc": pd.Timestamp("2019-01-01"),   # before AF3 2021 cutoff
               "2def": pd.Timestamp("2023-01-01"),   # after AF3 cutoff
               "3ghi": pd.Timestamp("2019-01-01")}   # after HistoFold 2018 cutoff
    out = classify(df, reldate)
    period = dict(zip(out.pdb_id, out.period))
    assert period["1abc"] == "Before"
    assert period["2def"] == "After"
    assert period["3ghi"] == "After"   # 2019 > 2018-04-30 HistoFold cutoff


def test_classify_unknown_when_no_date():
    df = pd.DataFrame({"method": ["Boltz-2"], "pdb_id": ["9zzz"], "rmsd_ca": [0.3]})
    out = classify(df, {})
    assert out.period.iloc[0] == "Unknown"


def test_summarize_ratio_and_significance():
    # 10 before (~high pLDDT), 10 after (clearly lower) -> ratio<1, small p
    import numpy as np
    rng = np.random.default_rng(0)
    before = pd.DataFrame({"method": "AlphaFold3", "pdb_id": [f"b{i}" for i in range(10)],
                           "period": "Before", "v": rng.normal(95, 1, 10)})
    after = pd.DataFrame({"method": "AlphaFold3", "pdb_id": [f"a{i}" for i in range(10)],
                          "period": "After", "v": rng.normal(88, 1, 10)})
    df = pd.concat([before, after], ignore_index=True)
    s = summarize(df, "v", "test").query("method=='AlphaFold3'").iloc[0]
    assert s["n_before"] == 10 and s["n_after"] == 10
    assert s["ratio_after_over_before"] < 1.0
    assert s["mwu_p"] < 0.05


def test_summarize_no_after_data_is_none():
    df = pd.DataFrame({"method": ["Boltz-2"] * 5, "pdb_id": [f"b{i}" for i in range(5)],
                       "period": "Before", "v": [95.0, 96, 97, 98, 99]})
    s = summarize(df, "v", "test").query("method=='Boltz-2'").iloc[0]
    assert s["n_after"] == 0
    assert s["mwu_p"] is None and s["ratio_after_over_before"] is None
