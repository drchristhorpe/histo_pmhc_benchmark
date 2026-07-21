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
@click.option("--no-figure", is_flag=True, help="Skip the summary figure.")
@click.option("--quiet", is_flag=True, help="Suppress progress output.")
def run(dataset: Path, out_dir: Path, aligned_root: Path | None,
        no_align: bool, no_figure: bool, quiet: bool) -> None:
    """Run the full pipeline on DATASET (a dir with predictions/, ground_truth/aligned/, metadata/)."""
    console = Console(quiet=quiet)
    console.print(f"[bold]pMHC benchmark[/bold] — dataset: {dataset}")
    manifest = run_pipeline(
        base=str(dataset), out_dir=str(out_dir),
        aligned_root=str(aligned_root) if aligned_root else None,
        do_align=not no_align, do_figure=not no_figure,
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


if __name__ == "__main__":
    main()
