# Open Cohort Factory

Open Cohort Factory builds reproducible disease cohorts and explicitly characterized reference
populations from public biomedical data.

Its central premise is that **normal is not a universal sample property**. Adjacent non-tumor
tissue, postmortem reference tissue, organ-donor tissue, and healthy-volunteer samples come from
different populations and collection processes. The software keeps those contexts visible in the
data model, provenance, and final report.

## Current capabilities

- Validate a human-readable YAML cohort specification.
- Query public case and sample metadata from the NCI Genomic Data Commons.
- Apply project, sample-type, age, sex, and primary-site filters.
- Normalize records into a source-independent schema.
- Preserve donor/sample nesting and tissue-acquisition context.
- Record the exact source query, retrieval time, and citation.
- Generate JSON manifests, summary statistics, and a self-contained HTML report.
- Materialize GTEx RNA-seq sample metadata, including age bracket, sex, RIN, ischemic time,
  Hardy scale, autolysis, and pathology annotations.
- Compare disease and reference populations at the donor level and recommend analysis safeguards.
- Materialize Tabula Sapiens metadata as donor × cell-type aggregates through a pinned CELLxGENE
  Census release.
- Report both cell counts and donor coverage for every cell type.
- Declare additional CELLxGENE reference panels in the same specification.

Single-cell support is optional because its local dependency stack is substantially larger than the
GDC and GTEx metadata clients.

## First workflow

```bash
uv sync --extra dev --extra single-cell
uv run cohort-factory validate examples/luad_older_adults.yaml
uv run cohort-factory build examples/luad_older_adults.yaml --output outputs/luad
```

Continue with the [installation guide](installation.md) or work through the
[lung adenocarcinoma vignette](vignettes/luad-older-adults.md).
