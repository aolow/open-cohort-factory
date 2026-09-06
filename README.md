# Open Cohort Factory

Open Cohort Factory builds reproducible disease cohorts and explicitly characterized reference
populations from public biomedical data. It treats “normal” as a claim that requires evidence,
not as a universal sample label.

The initial release retrieves public case and sample metadata from the NCI Genomic Data Commons,
normalizes it into a source-independent model, preserves tissue-acquisition context, and produces
an auditable HTML report. The schema already represents GTEx and Tabula Sapiens reference panels;
their data connectors are the next implementation milestone.

## Why this exists

Cancer cohorts and reference atlases differ in age, health history, tissue acquisition, processing,
assay, and metadata completeness. A postmortem GTEx sample, an organ-donor Tabula Sapiens sample,
and adjacent non-tumor tissue from a person with cancer answer related but different questions.
This project keeps those distinctions visible.

## Quick start

Requirements: Python 3.11 or newer and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --all-extras --dev
uv run cohort-factory validate examples/luad_older_adults.yaml
uv run cohort-factory build examples/luad_older_adults.yaml --output outputs/luad
open outputs/luad/report.html
```

## Outputs

- `manifest.json`: validated specification, normalized sample metadata, source query, and provenance
- `summary.json`: donor-aware counts and metadata-completeness summary
- `report.html`: a self-contained human-readable report with interpretation guardrails

## Reference contexts

The canonical model uses precise acquisition labels:

- `adjacent_non_tumor`
- `postmortem_reference`
- `organ_donor`
- `surgical_non_diseased`
- `healthy_volunteer`
- `disease_control`
- `unknown_or_incompletely_characterized`

Missing information is represented as unknown; it is never interpreted as absence of disease,
comorbidity, or exposure.

## Scientific guardrails

- Adjacent non-tumor tissue is not labeled healthy.
- Cross-study expression matrices are not treated as unconfounded differential-expression inputs.
- Donors, rather than cells, are the independent unit for inferential single-cell comparisons.
- Results remain stratified by source unless an explicit, documented harmonization method is used.
- Every build records its source query, retrieval time, citations, and known metadata gaps.

## Roadmap

1. GDC cohort metadata, provenance, validation, and reporting
2. GTEx bulk-tissue reference connector and age/sex comparability audit
3. Tabula Sapiens donor-aware cell-type reference connector
4. Balance diagnostics, matching or weighting, and sensitivity analyses across reference panels
5. Assay-aware expression summaries and an interactive cohort explorer

## Data use

The software is MIT licensed. Downloaded data remain governed by their original source terms and
are not committed to this repository. Users are responsible for reviewing the applicable data-use
policies and citations for each generated cohort.

