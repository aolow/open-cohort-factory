"""Command-line interface."""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from .audit import methodological_warnings, summarize
from .config import load_spec
from .models import BuildResult
from .report import write_outputs
from .sources.gdc import GDCClient
from .sources.gtex import GTExClient

app = typer.Typer(no_args_is_help=True, help="Build reproducible public-data cohorts.")
console = Console()


@app.command()
def validate(spec: Annotated[Path, typer.Argument(exists=True, readable=True)]) -> None:
    """Validate a cohort specification without downloading data."""
    project = load_spec(spec)
    console.print(f"[green]Valid[/green]: {project.name}")
    table = Table("Reference panel", "Source", "Context")
    for panel in project.reference_panels:
        table.add_row(panel.name, panel.source.value, panel.context.value)
    console.print(table)


@app.command()
def build(
    spec: Annotated[Path, typer.Argument(exists=True, readable=True)],
    output: Annotated[Path, typer.Option("--output", "-o")] = Path("outputs/latest"),
) -> None:
    """Materialize the disease cohort and create an auditable report."""
    project = load_spec(spec)
    with console.status("Retrieving public GDC metadata..."):
        samples, disease_provenance = GDCClient().fetch(project.disease_cohort)
    provenance = [disease_provenance]
    for panel in project.reference_panels:
        if panel.source.value != "gtex":
            continue
        with console.status(f"Retrieving public GTEx metadata for {panel.name}..."):
            reference_samples, reference_provenance = GTExClient().fetch(panel)
        samples.extend(reference_samples)
        provenance.append(reference_provenance)
    warnings = methodological_warnings(project, samples)
    result = BuildResult(
        specification=project, samples=samples, provenance=provenance, warnings=warnings
    )
    summary = summarize(samples, project)
    write_outputs(result, summary, output)
    console.print(
        f"[green]Built[/green] {summary['sample_count']} samples from "
        f"{summary['donor_count']} donors in {output}"
    )


if __name__ == "__main__":
    app()
