"""Click CLI for histo-pmhc-benchmark."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from histo_pmhc_benchmark.dataset import build_dataset
from histo_pmhc_benchmark.pipeline import run_pipeline
from histo_pmhc_benchmark import stages


def _progress(console: Console):
    def emit(msg: str):
        console.print(f"  [dim]·[/dim] {msg}")
    return emit


@click.group()
@click.version_option(package_name="histo-pmhc-benchmark")
def main() -> None:
    """Benchmark pMHC class I structure predictions against ground truth."""


@main.command()
@click.argument("dataset", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("-o", "--out-dir", type=click.Path(file_okay=False, path_type=Path),
              default="pmhc_benchmark_out", help="Output directory for CSVs + figure.")
@click.option("--aligned-root", type=click.Path(file_okay=False, path_type=Path), default=None,
              help="Where aligned PDBs live/are written (default: <out-dir>/aligned).")
@click.option("--no-align", is_flag=True, help="Reuse existing aligned PDBs (skip alignment).")
@click.option("--no-figure", is_flag=True, help="Skip the summary + box-plot figures.")
@click.option("--no-cutoff", is_flag=True, help="Skip the before/after training-cutoff analysis.")
@click.option("--offline", is_flag=True, help="Don't fetch missing release dates from RCSB/PDBe.")
@click.option("--quiet", is_flag=True, help="Suppress progress output.")
def run(dataset: Path, out_dir: Path, aligned_root: Path | None,
        no_align: bool, no_figure: bool, no_cutoff: bool, offline: bool, quiet: bool) -> None:
    """Run the full pipeline on DATASET (a dir with predictions/, ground_truth/aligned/, metadata/)."""
    console = Console(quiet=quiet)
    console.print(f"[bold]pMHC benchmark[/bold] — dataset: {dataset}")
    manifest = run_pipeline(
        base=str(dataset), out_dir=str(out_dir),
        aligned_root=str(aligned_root) if aligned_root else None,
        do_align=not no_align, do_figure=not no_figure,
        do_cutoff=not no_cutoff, fetch_missing_dates=not offline,
        progress=None if quiet else _progress(console),
    )
    if not quiet:
        t = Table(title="Headline medians (Å)")
        for col in ("Method", "n", "Peptide Cα", "Iface mainchain", "Iface sidechain"):
            t.add_column(col)
        for m, s in manifest["headline"].items():
            t.add_row(m, str(s["n"]), f"{s['peptide_ca_median']}",
                      f"{s['interface_mainchain_median']}", f"{s['interface_sidechain_median']}")
        console.print(t)
    console.print(f"[green]done[/green] — {len(manifest['files'])} files in {out_dir}")


@main.command()
@click.argument("dataset", type=click.Path(exists=True, file_okay=False, path_type=Path))
def index(dataset: Path) -> None:
    """Summarize a DATASET: GT complexes and per-method prediction counts."""
    console = Console()
    ds = build_dataset(str(dataset))
    t = Table(title=f"{dataset.name}: {len(ds.gt_ids)} GT complexes")
    t.add_column("Method"); t.add_column("Predicted"); t.add_column("Missing vs GT")
    for meth in __import__("histo_pmhc_benchmark.config", fromlist=["METHODS"]).METHODS:
        pred = set(ds.method_folders.get(meth.key, {}))
        t.add_row(meth.label, str(len(pred)), str(len(ds.gt_ids - pred)))
    console.print(t)


@main.command()
@click.argument("dataset", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("-o", "--out-dir", type=click.Path(file_okay=False, path_type=Path), default=".")
def gaps(dataset: Path, out_dir: Path) -> None:
    """Write per-method missing-complex CSVs (pdb_code,locus,allele_slug,peptide_sequence,resolution)."""
    console = Console()
    out_dir.mkdir(parents=True, exist_ok=True)
    ds = build_dataset(str(dataset))
    for slug, df in stages.gap_analysis(ds).items():
        p = out_dir / f"missing_{slug}.csv"
        df.to_csv(p, index=False)
        console.print(f"  {slug:14s} {len(df):4d} missing -> {p}")


@main.command()
@click.argument("dataset", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("-o", "--out-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
              default="pmhc_benchmark_out", help="Dir holding the scored *_whole.csv tables (from `run`).")
@click.option("--offline", is_flag=True, help="Don't fetch missing release dates from RCSB/PDBe.")
def cutoff(dataset: Path, out_dir: Path, offline: bool) -> None:
    """Before/after training-cutoff analysis on already-scored tables in OUT-DIR.

    Reads the peptide/interface *_whole.csv tables a prior `run` wrote, classifies
    each structure by its method's PDB cutoff, and writes the cutoff CSVs + box plots.
    """
    import pandas as pd
    from histo_pmhc_benchmark import cutoff as cutoff_mod
    from histo_pmhc_benchmark.figures import confidence_cutoff_boxplots

    console = Console()
    need = {"peptide_rmsd_whole": "pmhc_peptide_rmsd_whole.csv",
            "interface_rmsd_whole": "pmhc_interface_rmsd_whole.csv",
            "peptide_confidence_whole": "pmhc_peptide_confidence_whole.csv",
            "interface_confidence_whole": "pmhc_interface_confidence_whole.csv"}
    optional = {"peptide_confidence_per_residue": "pmhc_peptide_confidence_per_residue.csv",
                "interface_confidence_per_position": "pmhc_interface_confidence_per_position.csv"}
    tables = {}
    for k, fn in {**need, **optional}.items():
        p = out_dir / fn
        if p.exists():
            tables[k] = pd.read_csv(p)
        elif k in need:
            raise click.ClickException(f"missing required table {p} — run `histo-pmhc-benchmark run` first")

    res = cutoff_mod.run_cutoff_analysis(str(dataset), tables, fetch_missing=not offline,
                                         progress=lambda m: console.print(f"  [dim]·[/dim] {m}"))
    res["classification"].to_csv(out_dir / "pmhc_cutoff_classification.csv", index=False)
    res["accuracy_summary"].to_csv(out_dir / "pmhc_cutoff_accuracy_summary.csv", index=False)
    res["confidence_summary"].to_csv(out_dir / "pmhc_cutoff_confidence_summary.csv", index=False)
    confidence_cutoff_boxplots(res["peptide_confidence_classified"],
                               res["interface_confidence_classified"],
                               str(out_dir / "pmhc_confidence_cutoff_boxplots.png"))

    t = Table(title="Confidence (pLDDT) before → after cutoff (median)")
    for c in ("Metric", "Method", "n before", "n after", "before", "after", "MWU p"):
        t.add_column(c)
    for _, r in res["confidence_summary"].iterrows():
        t.add_row(r["metric"], r["method"], str(r["n_before"]), str(r["n_after"]),
                  str(r["median_before"]), str(r["median_after"]), str(r["mwu_p"]))
    console.print(t)
    console.print(f"[green]done[/green] — cutoff CSVs + box plots in {out_dir}")


if __name__ == "__main__":
    main()
