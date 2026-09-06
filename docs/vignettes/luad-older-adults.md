# Older adults with lung adenocarcinoma

This vignette builds a public TCGA lung adenocarcinoma cohort restricted to adults aged 60 or older,
materializes bulk postmortem lung metadata from GTEx, and declares a future cell-type-level
organ-donor lung reference from Tabula Sapiens.

The purpose is not yet to compare expression values. It is to determine whether the disease and
reference populations support that analysis and which adjustments or sensitivity analyses would
be required first.

## Inspect the specification

The complete configuration is in `examples/luad_older_adults.yaml`.

```yaml
disease_cohort:
  source: gdc
  projects:
    - TCGA-LUAD
  age:
    minimum: 60
  sample_types:
    - Primary Tumor

reference_panels:
  - name: GTEx lung older-adult reference
    source: gtex
    tissue: Lung
    context: postmortem_reference
    age:
      minimum: 60
      maximum: 70
    resolution: bulk

  - name: Tabula Sapiens lung cellular reference
    source: tabula_sapiens
    tissue: lung
    context: organ_donor
    age:
      minimum: 60
    resolution: cell_type
```

The age constraints are source-specific. They do not imply that the GDC disease cohort and either
reference panel are otherwise exchangeable.

## Validate before downloading

```bash
uv run cohort-factory validate examples/luad_older_adults.yaml
```

Validation checks enum values, required fields, and constraints such as the ordering of age bounds.
It performs no network request.

## Build the disease cohort

```bash
uv run cohort-factory build \
  examples/luad_older_adults.yaml \
  --output outputs/luad
```

The command queries public GDC metadata and writes:

```text
outputs/luad/
├── manifest.json
├── report.html
└── summary.json
```

The exact number of records may change as the source evolves. The report separates donor count from
sample count because a donor may contribute more than one qualifying sample.

## Interpret the current report

The report answers six initial questions:

1. How many qualifying donors and samples were found?
2. Are age and sex metadata present for the selected records?
3. Which reference panels are intended, and how was their tissue acquired?
4. Which methodological warnings must be resolved before quantitative comparison?
5. Which disease-cohort age strata are missing from the reference panel?
6. How different are the donor-level sex distributions and recorded specimen findings?

The GTEx panel appears as `materialized`; Tabula Sapiens remains
`declared_not_materialized`. Neither status implies that expression matrices have been harmonized.

## Follow the analysis decision

In a verified run on September 6, 2026, the configuration produced:

| Population | Donors | Age distribution | Female proportion |
| --- | ---: | --- | ---: |
| TCGA-LUAD disease cohort | 364 | 60–69: 46.7%; 70–79: 44.5%; 80–89: 8.8% | 52.7% |
| GTEx lung reference | 220 | 60–69: 90.0%; 70–79: 10.0% | 31.8% |

Counts may change as public sources evolve. The important result is the decision: the GTEx panel is
`context_only_for_full_cohort` because it does not cover the disease cohort's 80–89 age stratum and
its sex distribution is materially different.

The report therefore recommends:

1. Keep TCGA and GTEx summaries source-specific.
2. Restrict the primary matched analysis to covered ages or report older patients separately.
3. Stratify or weight descriptive comparisons by sex.
4. Examine sensitivity to GTEx pathology annotations.
5. Do not perform cross-study differential expression without assay-aware harmonization.

This is the intended workflow: cohort construction leads to an explicit analytical decision, not
automatically to a statistical test.

Generate a concise decision memo from any completed build:

```bash
uv run python examples/inspect_comparability.py outputs/luad/summary.json
```

This example consumes the machine-readable comparability result rather than reimplementing its
logic in a notebook. It can serve as the handoff between cohort construction and a downstream
analysis plan, workflow manager, or review document.

## Inspect specimen context

The GTEx connector preserves affirmative specimen-pathology categories and collection-related
variables. In the verified lung run, frequently recorded findings included congestion, emphysema,
fibrosis, edema, and hemorrhage. These observations show why `postmortem_reference` is more precise
than an unqualified `normal` label.

The manifest also retains RIN, ischemic time, Hardy death scale, autolysis score, tissue ontology,
and free-text pathology notes when supplied by the public API. Missing annotations remain unknown.

## What this prevents

The workflow prevents several quiet analytical errors:

- Calling adjacent tissue healthy.
- Treating absent comorbidity metadata as evidence that no comorbidity existed.
- Counting multiple samples from one person as independent donors.
- Presenting independently processed disease and reference matrices as an unconfounded contrast.
- Losing the exact query and retrieval time used to construct a cohort.

## Next extension

The next analytical extension will add an explicit matched-subcohort export rather than silently
discarding unmatched donors. The Tabula Sapiens adapter can then add donor-aware cell-type summaries
while keeping organ-donor context and study effects visible.
