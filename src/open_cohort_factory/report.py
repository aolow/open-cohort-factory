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
    :root { color-scheme: light; --ink:#14213d; --muted:#5f6b7a; --line:#dbe3ea; --bg:#f6f8fb; --accent:#137c8b; --tumor:#c84b61; --gtex:#2878b5; --adjacent:#8b5fbf; --cell:#2a9d78; }
    body { margin:0; background:var(--bg); color:var(--ink); font:16px/1.55 system-ui,-apple-system,sans-serif; }
    main { max-width:1040px; margin:0 auto; padding:48px 24px 80px; }
    h1 { font-size:2.2rem; margin:0 0 8px; } h2 { margin-top:40px; }
    .lede { color:var(--muted); max-width:760px; }
    .grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:16px; margin:28px 0; }
    .card { background:white; border:1px solid var(--line); border-radius:12px; padding:18px; }
    .metric { font-size:2rem; font-weight:700; color:var(--accent); }
    .label { color:var(--muted); font-size:.9rem; }
    table { width:100%; border-collapse:collapse; background:white; }
    th,td { border:1px solid var(--line); padding:10px 12px; text-align:left; vertical-align:top; overflow-wrap:anywhere; }
    th { background:#e8f1f3; } code { background:#e9edf2; padding:2px 5px; border-radius:4px; }
    .warning { border-left:4px solid #d97706; padding:8px 14px; margin:10px 0; background:#fffaf0; }
    .flow { display:grid; grid-template-columns:1fr auto 1fr auto minmax(240px,1.35fr); align-items:center; gap:12px; margin:28px 0; }
    .flow-node { background:white; border:1px solid var(--line); border-top:5px solid var(--accent); border-radius:12px; padding:16px; min-height:78px; }
    .flow-node strong { display:block; font-size:1.55rem; } .arrow { color:var(--muted); font-size:1.5rem; }
    .branches { display:grid; gap:9px; }.branch { background:white; border:1px solid var(--line); border-left:5px solid var(--gtex); border-radius:8px; padding:9px 12px; }
    .branch:nth-child(2) { border-left-color:var(--adjacent); }
    .context-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(215px,1fr)); gap:12px; margin:18px 0; }
    .context { background:white; border:1px solid var(--line); border-radius:12px; padding:15px; }
    .context .kind { display:inline-block; color:white; background:var(--gtex); border-radius:999px; padding:3px 9px; font-size:.78rem; font-weight:700; }
    .context:nth-child(2) .kind { background:var(--adjacent); }.context:nth-child(3) .kind { background:var(--cell); }
    .sensitivity { background:white; border:1px solid var(--line); border-radius:12px; padding:20px; margin:18px 0 28px; }
    .axis { display:flex; justify-content:space-between; color:var(--muted); font-size:.8rem; margin:0 0 12px 170px; }
    .gene-plot { display:grid; grid-template-columns:155px 1fr; gap:14px; margin:18px 0; align-items:start; }
    .gene-label { font-weight:750; padding-top:3px; }.estimate-list { border-left:1px solid var(--line); border-right:1px solid var(--line); }
    .estimate { display:grid; grid-template-columns:minmax(190px,1.2fr) minmax(180px,2fr) 62px; gap:10px; align-items:center; min-height:32px; font-size:.82rem; }
    .track { height:18px; position:relative; background:linear-gradient(to right,transparent 49.7%,var(--ink) 49.7%,var(--ink) 50.3%,transparent 50.3%); }
    .dot { position:absolute; top:4px; width:10px; height:10px; border-radius:50%; background:var(--gtex); transform:translateX(-50%); }
    .estimate:nth-child(2) .dot { background:#4d9acb; }.estimate:nth-child(3) .dot { background:var(--adjacent); }.estimate:nth-child(4) .dot { background:#b07bd3; }
    .value { text-align:right; font-variant-numeric:tabular-nums; font-weight:650; }
    .availability { min-width:105px; }.coverage-bar { height:8px; background:#edf1f5; border-radius:5px; overflow:hidden; margin-top:4px; }
    .coverage-fill { display:block; height:100%; background:var(--accent); }.confounded { color:#a64b00; font-weight:700; }
    .covariate-note { display:grid; grid-template-columns:repeat(auto-fit,minmax(230px,1fr)); gap:10px; margin:14px 0; }
    .covariate-note > div { border-left:4px solid #d97706; background:#fffaf0; padding:10px 13px; }
    .ci { position:absolute; top:8px; height:2px; background:var(--ink); opacity:.55; }
    .fitness-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(285px,1fr)); gap:16px; margin:18px 0; }
    .fitness-panel { background:white; border:1px solid var(--line); border-radius:12px; padding:17px; }
    .dimension { display:grid; grid-template-columns:115px 1fr; gap:9px; padding:9px 0; border-top:1px solid var(--line); }
    .status { font-size:.72rem; font-weight:800; text-transform:uppercase; letter-spacing:.04em; }
    .status-supported { color:#18724e; }.status-caution { color:#a05a00; }.status-limited { color:#b4233b; }.status-not_evaluable { color:var(--muted); }
    .use-columns { display:grid; grid-template-columns:1fr 1fr; gap:12px; font-size:.86rem; }.use-columns ul { padding-left:18px; }
    .safety-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(290px,1fr)); gap:16px; margin:18px 0; }
    .safety-card { background:white; border:1px solid var(--line); border-radius:12px; padding:17px; }
    .tissue-row { display:grid; grid-template-columns:105px 1fr 42px; gap:8px; align-items:center; margin:7px 0; font-size:.8rem; }
    .tissue-track { background:#edf1f5; height:9px; border-radius:6px; overflow:hidden; }.tissue-fill { display:block; height:100%; background:var(--gtex); }
    @media (max-width:720px) { .flow { grid-template-columns:1fr; }.arrow { transform:rotate(90deg); text-align:center; }.axis { margin-left:0; }.gene-plot { grid-template-columns:1fr; }.estimate { grid-template-columns:1fr; gap:2px; padding:5px 0; } }
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
  {% if summary.expression is defined %}
  <h2>Analysis flow</h2>
  <div class="flow" aria-label="Expression cohort attrition">
    <div class="flow-node"><strong>{{ summary.expression.coverage.requested_sample_count }}</strong><span class="label">eligible cohort identifiers requested</span></div>
    <div class="arrow">→</div>
    <div class="flow-node"><strong>{{ summary.expression.coverage.returned_sample_count }}</strong><span class="label">samples represented in Xena Toil</span></div>
    <div class="arrow">→</div>
    <div class="branches">
      {% for analysis in summary.expression.comparisons %}{% if analysis.complete_expression_pairs is defined %}<div class="branch"><strong>{{ analysis.complete_expression_pairs }} pairs</strong> · demographic match with expression</div>{% endif %}{% if analysis.within_donor_pairs is defined %}<div class="branch"><strong>{{ analysis.within_donor_pairs }} pairs</strong> · within-participant tumor/adjacent</div>{% endif %}{% endfor %}
    </div>
  </div>
  {% endif %}
  <h2>Reference panels</h2>
  <div class="context-grid">
  {% for panel in summary.reference_panels %}<div class="context"><span class="kind">{{ panel.context|replace("_", " ") }}</span><h3>{{ panel.name }}</h3><p>{{ panel.source }} · {{ panel.resolution }} · {{ panel.tissue }}</p></div>{% endfor %}
  </div>
  <table><thead><tr><th>Name</th><th>Source</th><th>Tissue</th><th>Context</th><th>Resolution</th><th>Status</th></tr></thead>
  <tbody>{% for panel in summary.reference_panels %}<tr><td>{{ panel.name }}</td><td>{{ panel.source }}</td><td>{{ panel.tissue }}</td><td><code>{{ panel.context }}</code></td><td>{{ panel.resolution }}</td><td>{{ panel.status }}</td></tr>{% endfor %}</tbody></table>
  <h2>Population comparison</h2>
  <table><thead><tr><th>Role</th><th>Population</th><th>Source</th><th>Records</th><th>Donors</th><th>Cells</th><th>Age brackets</th><th>Sex</th></tr></thead>
  <tbody>{% for population in summary.populations %}<tr><td>{{ population.role }}</td><td>{{ population.name }}</td><td>{{ population.source }}</td><td>{{ population.samples }}</td><td>{{ population.donors }}</td><td>{{ population.cells if population.cells is not none else "—" }}</td><td>{{ population.age_brackets }}</td><td>{{ population.sex }}</td></tr>{% endfor %}</tbody></table>
  <h2>Covariate landscape</h2>
  <p>{{ summary.covariate_landscape.interpretation }}</p>
  <table><thead><tr><th>Covariate</th>{% for population in summary.covariate_landscape.populations %}<th>{{ population }}</th>{% endfor %}<th>Assessment</th></tr></thead>
  <tbody>{% for field in summary.covariate_landscape.fields %}<tr><td><strong>{{ field.label }}</strong><br><small>{{ field.role }}</small></td>{% for item in field.coverage %}<td class="availability"><strong>{{ item.percent }}%</strong> known<div class="coverage-bar"><span class="coverage-fill" style="width:{{ item.percent }}%"></span></div>{% if item.detail.median is defined %}<small>median {{ item.detail.median }}; {{ item.detail.minimum }}–{{ item.detail.maximum }}</small>{% elif item.detail.top_values is defined %}<small>{% for value in item.detail.top_values %}{{ value.value }} ({{ value.donors }}){% if not loop.last %}; {% endif %}{% endfor %}</small>{% endif %}</td>{% endfor %}<td>{% if field.availability_status == "source_confounded" %}<span class="confounded">source-confounded availability</span>{% elif field.availability_status == "unavailable" %}unavailable in all panels{% elif field.distribution_shift %}<span class="confounded">distribution differs across populations</span>{% else %}no large observed shift{% endif %}</td></tr>{% endfor %}</tbody></table>
  <div class="covariate-note">{% for field in summary.covariate_landscape.fields %}{% if field.source_confounded_availability %}<div><strong>{{ field.label }}</strong><br>Availability differs sharply across populations; adjustment may select a source-specific subset.</div>{% elif field.distribution_shift %}<div><strong>{{ field.label }}</strong><br>Observed distributions differ across sufficiently represented populations (maximum distance {{ field.maximum_distribution_distance }}).</div>{% endif %}{% endfor %}</div>
  <h2>Reference fitness by intended use</h2>
  <p>No composite score is calculated: strength in one dimension cannot cancel a structural limitation in another.</p>
  <div class="fitness-grid">{% for panel in summary.reference_fitness %}<section class="fitness-panel"><h3>{{ panel.reference_panel }}</h3><p><code>{{ panel.context }}</code></p>
    {% for item in panel.dimensions %}<div class="dimension"><div><span class="status status-{{ item.status }}">{{ item.status|replace("_", " ") }}</span><br>{{ item.dimension }}</div><div>{{ item.evidence }}</div></div>{% endfor %}
    <div class="use-columns"><div><strong>Appropriate uses</strong><ul>{% for use in panel.recommended_uses %}<li>{{ use }}</li>{% endfor %}</ul></div><div><strong>Not supported</strong><ul>{% for use in panel.not_supported %}<li>{{ use }}</li>{% endfor %}</ul></div></div>
  </section>{% endfor %}</div>
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
    {% if match.method == "within_donor" %}<p><strong>{{ match.matched_pairs }}</strong> tumor/reference pairs linked within the same participant.</p>
    {% else %}<p><strong>{{ match.matched_pairs }}</strong> donor pairs matched exactly on {{ match.variables|join(", ") }} using seed {{ match.seed }}.</p>{% endif %}
    <table><thead><tr>{% for variable in match.variables %}<th>{{ variable }}</th>{% endfor %}<th>Matched pairs</th></tr></thead>
    <tbody>{% for row in match.matched_strata %}<tr>{% for variable in match.variables %}<td>{{ row.stratum[variable] }}</td>{% endfor %}<td>{{ row.matched_pairs }}</td></tr>{% endfor %}</tbody></table>
    <p>{{ match.excluded_donors|length }} donor records were not selected. Reasons are retained in <code>matching.json</code>.</p>
  </div>
  {% else %}<p>Matching is disabled or no reference panel has been materialized.</p>{% endfor %}
  <h2>Expression comparison</h2>
  {% if summary.expression is defined %}
  <p>Values are cohort-restricted UCSC Xena Toil measurements on the {{ summary.expression.scale }} scale. Differences are descriptive and do not imply causal or universally healthy reference populations.</p>
  <div class="sensitivity" role="img" aria-label="Median tumor minus reference expression differences across reference definitions">
    <h3>How the reference definition changes the signal</h3>
    <div class="axis"><span>{{ summary.expression.sensitivity_chart.minimum }}</span><span>tumor − reference median difference</span><span>+{{ summary.expression.sensitivity_chart.maximum }}</span></div>
    {% for gene in summary.expression.sensitivity_chart.genes %}<div class="gene-plot"><div class="gene-label">{{ gene.gene }}</div><div class="estimate-list">
      {% for estimate in gene.estimates %}<div class="estimate"><span>{{ estimate.comparison }}</span><div class="track"><span class="ci" style="left:{{ estimate.ci_left }}%;width:{{ estimate.ci_width }}%"></span><span class="dot" style="left:{{ estimate.position }}%"></span></div><span class="value">{{ estimate.value }}<br><small>{{ estimate.ci_lower }}, {{ estimate.ci_upper }}</small></span></div>{% endfor %}
    </div></div>{% endfor %}
  </div>
  {% for analysis in summary.expression.comparisons %}
  <div class="card"><h3>{{ analysis.name }}</h3>
    {% if analysis.complete_expression_pairs is defined %}<p><strong>{{ analysis.complete_expression_pairs }}</strong> matched pairs had expression available for both donors.</p>{% endif %}
    {% if analysis.within_donor_pairs is defined %}<p><strong>{{ analysis.within_donor_pairs }}</strong> participants contributed both tumor and adjacent non-tumor expression.</p>{% endif %}
    <table><thead><tr><th>Gene</th><th>Disease donors</th><th>Reference donors</th><th>Disease median (IQR)</th><th>Reference median (IQR)</th><th>Median difference (95% bootstrap CI)</th></tr></thead>
    <tbody>{% for row in analysis.genes %}<tr><td>{{ row.gene }}</td><td>{{ row.disease_donors }}</td><td>{{ row.reference_donors }}</td><td>{{ row.disease_median }} ({{ row.disease_q1 }}–{{ row.disease_q3 }})</td><td>{{ row.reference_median }} ({{ row.reference_q1 }}–{{ row.reference_q3 }})</td><td>{{ row.median_difference }} ({{ row.ci_lower }}, {{ row.ci_upper }})</td></tr>{% else %}<tr><td colspan="6">No comparable Xena measurements were available.</td></tr>{% endfor %}</tbody></table>
  </div>
  {% endfor %}
  {% if summary.expression.pan_tissue is defined %}
  <h2>Pan-tissue normal expression</h2>
  <p>GTEx tissues are ranked independently for each gene. Bars show median {{ summary.expression.pan_tissue.scale }} relative to that gene's highest-expression tissue; prevalence uses TPM ≥ {{ summary.expression.pan_tissue.tpm_threshold }}.</p>
  <div class="safety-grid">{% for gene in summary.expression.pan_tissue.genes %}<section class="safety-card"><h3>{{ gene.gene }}</h3><p><strong>{{ gene.tissues_with_majority_above_threshold }}</strong> tissues have at least half of donors above the configured threshold.</p>
    {% for tissue in gene.top_tissues %}<div class="tissue-row"><span>{{ tissue.tissue }}</span><div class="tissue-track"><span class="tissue-fill" style="width:{{ tissue.chart_percent }}%"></span></div><span>{{ tissue.median }}</span></div>{% endfor %}
  </section>{% endfor %}</div>
  {% endif %}
  {% if summary.expression.cell_type_attribution is defined %}
  <h2>Which normal cells carry the signal?</h2>
  <p>Tabula Sapiens raw counts are summarized within donor × cell type. Detection is descriptive—not evidence that a transcript is surface-accessible or that targeting it causes toxicity.</p>
  {% for panel in summary.expression.cell_type_attribution %}
  <details><summary><strong>{{ panel.reference_panel }}</strong> · donor-aware attribution</summary>
  {% for gene in panel.genes %}<h3>{{ gene.gene }}</h3><div class="table-wrap"><table><thead><tr><th>Cell type</th><th>Donors</th><th>Cells</th><th>Median donor detection</th><th>Median mean log1p count</th></tr></thead>
  <tbody>{% for row in gene.top_cell_types %}<tr><td>{{ row.cell_type }}</td><td>{{ row.donors }}</td><td>{{ row.cells }}</td><td>{{ (row.median_fraction_detected * 100)|round(1) }}%</td><td>{{ row.median_mean_log1p_raw_count }}</td></tr>{% endfor %}</tbody></table></div>{% endfor %}
  </details>
  {% endfor %}
  {% endif %}
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
  <details><summary><strong>Quality notes, metadata, and provenance</strong></summary>
  <h2>Interpretation guardrails</h2>
  {% for warning in result.warnings %}<div class="warning">{{ warning }}</div>{% endfor %}
  <h2>Metadata completeness</h2>
  <table><thead><tr><th>Field</th><th>Missing samples</th></tr></thead>
  <tbody>{% for field,count in summary.missing_fields.items() %}<tr><td>{{ field }}</td><td>{{ count }}</td></tr>{% else %}<tr><td colspan="2">No required metadata fields are missing.</td></tr>{% endfor %}</tbody></table>
  <h2>Provenance</h2>
  {% for item in result.provenance %}<p><strong>{{ item.source.value }}</strong> — retrieved {{ item.retrieved_at }} from <a href="{{ item.citation_url }}">source documentation</a>.</p>{% endfor %}
  </details>
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
            fieldnames = list(dict.fromkeys(key for row in expression_rows for key in row))
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(expression_rows)
    attribution_rows = []
    for panel in summary.get("expression", {}).get("cell_type_attribution", []):
        for gene in panel["genes"]:
            for row in gene["cell_types"]:
                attribution_rows.append(
                    {
                        "reference_panel": panel["reference_panel"],
                        "gene": gene["gene"],
                        "top_cell_type": gene["top_cell_type"],
                        "leave_one_donor_out_top_agreement": gene[
                            "leave_one_donor_out_top_agreement"
                        ],
                        **row,
                    }
                )
    if attribution_rows:
        with (output_dir / "cell_type_attribution.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(attribution_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(attribution_rows)
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
