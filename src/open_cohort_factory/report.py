"""Self-contained HTML reporting."""

from __future__ import annotations

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
    <div class="card"><div class="metric">{{ summary.sample_count }}</div><div class="label">samples</div></div>
    <div class="card"><div class="metric">{{ summary.donor_count }}</div><div class="label">donors</div></div>
    <div class="card"><div class="metric">{{ summary.age.mean if summary.age.mean is not none else "—" }}</div><div class="label">mean age with metadata</div></div>
    <div class="card"><div class="metric">{{ summary.reference_panels|length }}</div><div class="label">declared reference panels</div></div>
  </div>
  <h2>Reference panels</h2>
  <table><thead><tr><th>Name</th><th>Source</th><th>Tissue</th><th>Context</th><th>Resolution</th><th>Status</th></tr></thead>
  <tbody>{% for panel in summary.reference_panels %}<tr><td>{{ panel.name }}</td><td>{{ panel.source }}</td><td>{{ panel.tissue }}</td><td><code>{{ panel.context }}</code></td><td>{{ panel.resolution }}</td><td>{{ panel.status }}</td></tr>{% endfor %}</tbody></table>
  <h2>Interpretation guardrails</h2>
  {% for warning in result.warnings %}<div class="warning">{{ warning }}</div>{% endfor %}
  <h2>Metadata completeness</h2>
  <table><thead><tr><th>Field</th><th>Missing samples</th></tr></thead>
  <tbody>{% for field,count in summary.missing_fields.items() %}<tr><td>{{ field }}</td><td>{{ count }}</td></tr>{% else %}<tr><td colspan="2">No required metadata fields are missing.</td></tr>{% endfor %}</tbody></table>
  <h2>Provenance</h2>
  {% for item in result.provenance %}<p><strong>{{ item.source.value }}</strong> — retrieved {{ item.retrieved_at }} from <a href="{{ item.citation_url }}">source documentation</a>.</p>{% endfor %}
</main></body></html>"""


def write_outputs(result: BuildResult, summary: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "manifest.json").write_text(
        json.dumps(result.model_dump(mode="json"), indent=2), encoding="utf-8"
    )
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    html = Environment(loader=BaseLoader(), autoescape=True).from_string(TEMPLATE).render(
        result=result, summary=summary
    )
    (output_dir / "report.html").write_text(html, encoding="utf-8")

