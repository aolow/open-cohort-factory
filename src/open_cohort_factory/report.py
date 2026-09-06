"""Self-contained HTML reporting."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from jinja2 import BaseLoader, Environment

from .models import BuildResult

TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{{ result.specification.name }} cohort report</title>
  <style>
    :root { color-scheme: light; --ink:#14213d; --muted:#5f6b7a; --line:#dbe3ea; --bg:#f6f8fb; --accent:#137c8b; }
    body { margin:0; background:var(--bg); color:var(--ink); font:16px/1.55 system-ui,-apple-system,sans-serif; }
    main { max-width:1040px; margin:0 auto; padding:48px 24px 80px; }
    h1 { font-size:2.2rem; margin:0 0 8px; } h2 { margin-top:40px; }
    .lede { color:var(--muted); max-width:760px; }
    .grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:16px; margin:28px 0; }
    .card { background:white; border:1px solid var(--line); border-radius:12px; padding:18px; }
    .metric { font-size:2rem; font-weight:700; color:var(--accent); }
    .label { color:var(--muted); font-size:.9rem; }
    table { width:100%; border-collapse:collapse; background:white; }
    th,td { border:1px solid var(--line); padding:10px 12px; text-align:left; }
    th { background:#e8f1f3; } code { background:#e9edf2; padding:2px 5px; border-radius:4px; }
    .warning { border-left:4px solid #d97706; padding:8px 14px; margin:10px 0; background:#fffaf0; }
  </style>
</head>
<body><main>
  <h1>{{ result.specification.name }}</h1>
  <p class="lede">{{ result.specification.description }}</p>
  <div class="grid">
    <div class="card"><div class="metric">{{ summary.sample_count }}</div><div class="label">normalized records</div></div>
    <div class="card"><div class="metric">{{ summary.donor_count }}</div><div class="label">donors</div></div>
    <div class="card"><div class="metric">{{ summary.age.mean if summary.age.mean is not none else "—" }}</div><div class="label">disease-cohort mean age</div></div>
    <div class="card"><div class="metric">{{ summary.reference_panels|length }}</div><div class="label">declared reference panels</div></div>
  </div>
  <h2>Reference panels</h2>
  <table><thead><tr><th>Name</th><th>Source</th><th>Tissue</th><th>Context</th><th>Resolution</th><th>Status</th></tr></thead>
  <tbody>{% for panel in summary.reference_panels %}<tr><td>{{ panel.name }}</td><td>{{ panel.source }}</td><td>{{ panel.tissue }}</td><td><code>{{ panel.context }}</code></td><td>{{ panel.resolution }}</td><td>{{ panel.status }}</td></tr>{% endfor %}</tbody></table>
  <h2>Population comparison</h2>
  <table><thead><tr><th>Role</th><th>Population</th><th>Source</th><th>Records</th><th>Donors</th><th>Cells</th><th>Age brackets</th><th>Sex</th></tr></thead>
  <tbody>{% for population in summary.populations %}<tr><td>{{ population.role }}</td><td>{{ population.name }}</td><td>{{ population.source }}</td><td>{{ population.samples }}</td><td>{{ population.donors }}</td><td>{{ population.cells if population.cells is not none else "—" }}</td><td>{{ population.age_brackets }}</td><td>{{ population.sex }}</td></tr>{% endfor %}</tbody></table>
  <h2>Comparability assessment</h2>
  {% for comparison in summary.comparability %}
  <div class="card">
    <h3>{{ comparison.reference_panel }}</h3>
    <p><strong>Decision:</strong> <code>{{ comparison.decision }}</code></p>
    <p><strong>Disease age brackets:</strong> {{ comparison.disease_age_brackets }}</p>
    <p><strong>Reference age brackets:</strong> {{ comparison.reference_age_brackets }}</p>
    <p><strong>Female proportion:</strong> disease {{ comparison.female_proportion.disease }}, reference {{ comparison.female_proportion.reference }}, absolute difference {{ comparison.female_proportion.absolute_difference }}</p>
    <p><strong>Recommended actions</strong></p>
    <ul>{% for action in comparison.recommended_actions %}<li>{{ action }}</li>{% endfor %}</ul>
    <p><strong>Most frequent recorded specimen findings</strong></p>
    <ul>{% for item in comparison.top_reference_pathology_categories %}<li>{{ item.category }}: {{ item.samples }} samples</li>{% else %}<li>No affirmative pathology categories recorded.</li>{% endfor %}</ul>
  </div>
  {% else %}<p>No reference panel has been materialized.</p>{% endfor %}
  <h2>Matched cohorts</h2>
  {% for match in summary.matching %}
  <div class="card">
    <h3>{{ match.reference_panel }}</h3>
    <p><strong>{{ match.matched_pairs }}</strong> donor pairs matched exactly on {{ match.variables|join(", ") }} using seed {{ match.seed }}.</p>
    <table><thead><tr>{% for variable in match.variables %}<th>{{ variable }}</th>{% endfor %}<th>Matched pairs</th></tr></thead>
    <tbody>{% for row in match.matched_strata %}<tr>{% for variable in match.variables %}<td>{{ row.stratum[variable] }}</td>{% endfor %}<td>{{ row.matched_pairs }}</td></tr>{% endfor %}</tbody></table>
    <p>{{ match.excluded_donors|length }} donor records were not selected. Reasons are retained in <code>matching.json</code>.</p>
  </div>
  {% else %}<p>Matching is disabled or no reference panel has been materialized.</p>{% endfor %}
  <h2>Expression comparison</h2>
  {% if summary.expression is defined %}
  <p>Values are cohort-restricted UCSC Xena Toil measurements on the {{ summary.expression.scale }} scale. Differences are descriptive and do not imply causal or universally healthy reference populations.</p>
  {% for analysis in summary.expression.comparisons %}
  <div class="card"><h3>{{ analysis.name }}</h3>
    {% if analysis.complete_expression_pairs is defined %}<p><strong>{{ analysis.complete_expression_pairs }}</strong> matched pairs had expression available for both donors.</p>{% endif %}
    <table><thead><tr><th>Gene</th><th>Disease samples</th><th>Reference samples</th><th>Disease median</th><th>Reference median</th><th>Median difference</th></tr></thead>
    <tbody>{% for row in analysis.genes %}<tr><td>{{ row.gene }}</td><td>{{ row.disease_samples }}</td><td>{{ row.reference_samples }}</td><td>{{ row.disease_median }}</td><td>{{ row.reference_median }}</td><td>{{ row.median_difference }}</td></tr>{% else %}<tr><td colspan="6">No comparable Xena measurements were available.</td></tr>{% endfor %}</tbody></table>
  </div>
  {% endfor %}
  {% else %}<p>No expression retrieval was requested.</p>{% endif %}
  <h2>Single-cell reference context</h2>
  {% for panel in summary.cell_type_summaries %}
  <div class="card">
    <h3>{{ panel.reference_panel }}</h3>
    <p>{{ panel.total_cells }} cells summarized across {{ panel.total_donors }} eligible donors. Cell counts describe sampling depth; donor coverage determines the population evidence.</p>
    <table><thead><tr><th>Cell type</th><th>Ontology</th><th>Cells</th><th>Donors represented</th></tr></thead>
    <tbody>{% for row in panel.cell_types[:20] %}<tr><td>{{ row.cell_type }}</td><td><code>{{ row.ontology_term_id }}</code></td><td>{{ row.cells }}</td><td>{{ row.donors }}</td></tr>{% endfor %}</tbody></table>
  </div>
  {% else %}<p>No single-cell reference panel has been materialized.</p>{% endfor %}
  <h2>Interpretation guardrails</h2>
  {% for warning in result.warnings %}<div class="warning">{{ warning }}</div>{% endfor %}
  <h2>Metadata completeness</h2>
  <table><thead><tr><th>Field</th><th>Missing samples</th></tr></thead>
  <tbody>{% for field,count in summary.missing_fields.items() %}<tr><td>{{ field }}</td><td>{{ count }}</td></tr>{% else %}<tr><td colspan="2">No required metadata fields are missing.</td></tr>{% endfor %}</tbody></table>
  <h2>Provenance</h2>
  {% for item in result.provenance %}<p><strong>{{ item.source.value }}</strong> — retrieved {{ item.retrieved_at }} from <a href="{{ item.citation_url }}">source documentation</a>.</p>{% endfor %}
</main></body></html>"""


def write_outputs(
    result: BuildResult,
    summary: dict[str, Any],
    output_dir: Path,
    expression_rows: list[dict[str, Any]] | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "manifest.json").write_text(
        json.dumps(result.model_dump(mode="json"), indent=2), encoding="utf-8"
    )
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    expression_rows = expression_rows or []
    if expression_rows:
        with (output_dir / "expression.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(expression_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(expression_rows)
    matching = summary.get("matching", [])
    (output_dir / "matching.json").write_text(json.dumps(matching, indent=2), encoding="utf-8")
    with (output_dir / "matched_donors.tsv").open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "reference_panel",
            "pair_id",
            "source",
            "cohort_role",
            "cohort_name",
            "donor_id",
            "stratum",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for match in matching:
            for row in match["selected_donors"]:
                writer.writerow(
                    {
                        "reference_panel": match["reference_panel"],
                        **row,
                        "stratum": json.dumps(row["stratum"], sort_keys=True),
                    }
                )
    html = (
        Environment(loader=BaseLoader(), autoescape=True)
        .from_string(TEMPLATE)
        .render(result=result, summary=summary)
    )
    (output_dir / "report.html").write_text(html, encoding="utf-8")
