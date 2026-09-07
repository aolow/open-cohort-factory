# Open Cohort Factory

Open Cohort Factory builds reproducible disease cohorts and explicitly characterized reference
populations from public biomedical data. It treats “normal” as a claim that requires evidence,
not as a universal sample label.

The current release retrieves public case and sample metadata from the NCI Genomic Data Commons and
GTEx, donor-aware cell metadata from Tabula Sapiens through the CELLxGENE Census, and optional
cohort-restricted expression slices from the UCSC Xena Toil recompute. It preserves
tissue-acquisition context and produces an auditable HTML report.

## Why this exists

Cancer cohorts and reference atlases differ in age, health history, tissue acquisition, processing,
assay, and metadata completeness. A postmortem GTEx sample, an organ-donor Tabula Sapiens sample,
and adjacent non-tumor tissue from a person with cancer answer related but different questions.
This project keeps those distinctions visible.

## Quick start

Requirements: Python 3.11 or newer and [uv](https://docs.astral.sh/uv/).

```bash
uv sync --extra dev
uv run cohort-factory validate examples/luad_older_adults.yaml
uv run cohort-factory build examples/luad_older_adults.yaml --output outputs/luad
open outputs/luad/report.html
```

To materialize single-cell reference panels such as Tabula Sapiens, install the optional stack with
`uv sync --extra dev --extra single-cell`.

For complete setup instructions, concepts, and a worked example, see the `docs/` directory. Build
the local documentation site with `uv run mkdocs serve`.

## Outputs

- `manifest.json`: validated specification, normalized sample metadata, source query, and provenance
- `summary.json`: donor-aware counts and metadata-completeness summary
- `report.html`: a self-contained human-readable report with interpretation guardrails
- `matching.json`: selected and excluded donors with strata and explicit reasons
- `matched_donors.tsv`: analysis-ready table of deterministic donor pairs
- `expression.tsv`: requested Xena measurements only; never the full expression atlas

The HTML report includes a cross-population covariate-availability matrix, source-confounding
warnings, expression distributions, and deterministic 95% bootstrap intervals. Matched and
within-participant analyses resample donor pairs together.

Each reference panel also receives a multi-dimensional fitness assessment covering population
support, demographic alignment, acquisition, design linkage, measurement overlap, and baseline
metadata. The tool deliberately does not collapse these dimensions into a single score.

Optional pan-tissue mode profiles requested genes across every GTEx primary tissue represented in
the Xena Toil compendium. It reports donor-level medians, dispersion, and prevalence above a declared
TPM threshold so tissue-specific enrichment is not mistaken for whole-body selectivity.

Optional cell-type attribution then queries the same genes in declared Tabula Sapiens panels and
summarizes raw counts within donor × cell type. This localizes a bulk-tissue signal without treating
individual cells as independent biological replicates or equating RNA detection with target safety.

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
- Exact matching is performed at the donor level and never silently drops unmatched donors.
- Values such as `unknown` and `not reported` count as missing, not as observed categories.

## Roadmap

1. GDC cohort metadata, provenance, validation, and reporting
2. GTEx bulk-tissue reference connector and donor-level age/sex comparability audit
3. Tabula Sapiens donor-aware cell-type reference connector
4. Matching and auditable sensitivity analyses across reference panels
5. Cohort-restricted Xena Toil expression summaries without downloading whole-atlas matrices
6. Sensitivity analysis across reference definitions and an interactive cohort explorer

## Data use

The software is MIT licensed. Downloaded data remain governed by their original source terms and
are not committed to this repository. Users are responsible for reviewing the applicable data-use
policies and citations for each generated cohort.
