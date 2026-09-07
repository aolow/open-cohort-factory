"""Command-line interface."""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from .audit import methodological_warnings, summarize
from .config import load_spec
from .expression import summarize_expression, summarize_pan_tissue
from .matching import exact_match
from .models import BuildResult
from .reference_fitness import assess_reference_fitness
from .report import write_outputs
from .sources.gdc import GDCClient
from .sources.gtex import GTExClient
from .sources.tabula_sapiens import TabulaSapiensClient
from .sources.xena import XenaClient

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
        if panel.source.value == "gdc":
            with console.status(f"Retrieving public GDC metadata for {panel.name}..."):
                reference_samples, reference_provenance = GDCClient().fetch_reference(
                    panel, project.disease_cohort
                )
        elif panel.source.value == "gtex":
            with console.status(f"Retrieving public GTEx metadata for {panel.name}..."):
                reference_samples, reference_provenance = GTExClient().fetch(panel)
        elif panel.source.value == "tabula_sapiens":
            with console.status(f"Retrieving Tabula Sapiens cell metadata for {panel.name}..."):
                reference_samples, reference_provenance = TabulaSapiensClient().fetch(panel)
        else:
            continue
        samples.extend(reference_samples)
        provenance.append(reference_provenance)
    matching = exact_match(samples, project.matching)
    expression_rows: list[dict[str, object]] = []
    expression_summary: dict[str, object] | None = None
    if project.expression is not None:
        with console.status("Retrieving a cohort-restricted UCSC Xena expression slice..."):
            expression_rows, expression_provenance = XenaClient().fetch(
                samples, project.expression
            )
        provenance.append(expression_provenance)
        expression_summary = summarize_expression(expression_rows, matching)
        expression_summary["coverage"] = {
            key: expression_provenance.query[key]
            for key in (
                "requested_sample_count",
                "returned_sample_count",
                "returned_samples_by_source",
            )
        }
        if project.expression.gtex_pan_tissue:
            with console.status("Profiling selected genes across GTEx tissues in Xena..."):
                atlas_rows, atlas_provenance = XenaClient().fetch_gtex_atlas(project.expression)
            expression_rows.extend(atlas_rows)
            expression_summary["pan_tissue"] = summarize_pan_tissue(
                atlas_rows, project.expression.normal_tissue_tpm_threshold
            )
            expression_provenance.query["gtex_pan_tissue"] = atlas_provenance
    warnings = methodological_warnings(project, samples)
    result = BuildResult(
        specification=project, samples=samples, provenance=provenance, warnings=warnings
    )
    summary = summarize(samples, project)
    summary["matching"] = matching
    if expression_summary is not None:
        summary["expression"] = expression_summary
    summary["reference_fitness"] = assess_reference_fitness(
        samples, summary, matching, expression_summary
    )
    write_outputs(result, summary, output, expression_rows)
    console.print(
        f"[green]Built[/green] {summary['sample_count']} normalized records from "
        f"{summary['donor_count']} donors in {output}"
    )


if __name__ == "__main__":
    app()
