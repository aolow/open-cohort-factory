# Development

## Quality checks

Run the same checks expected for every contribution:

```bash
uv run ruff check .
uv run mypy src
uv run pytest
uv run mkdocs build --strict
```

## Repository layout

```text
src/open_cohort_factory/
├── audit.py       # donor-aware summaries and methodological warnings
├── comparability.py # population-overlap decisions and recommendations
├── cli.py         # command-line interface
├── config.py      # YAML loading and validation
├── models.py      # canonical data contracts
├── report.py      # self-contained HTML output
└── sources/
    ├── gdc.py     # public GDC adapter
    └── gtex.py    # public GTEx Portal API adapter
```

## Adding a data source

A connector should translate source records into canonical sample and provenance models while
retaining source-specific acquisition context. It must include mocked unit tests, a public-data
citation, an explicit data-use note, and at least one end-to-end example specification.

Do not commit downloaded biomedical datasets or generated cohort outputs.
