# Older adults with lung adenocarcinoma

This vignette builds a public TCGA lung adenocarcinoma cohort restricted to adults aged 60 or older.
It also declares two future reference panels: bulk postmortem lung from GTEx and a cell-type-level
organ-donor lung reference from Tabula Sapiens.

The purpose is not yet to compare expression values. It is to establish an auditable disease cohort
and make the intended reference populations—and their limitations—explicit before analysis begins.

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

The report answers four initial questions:

1. How many qualifying donors and samples were found?
2. Are age and sex metadata present for the selected records?
3. Which reference panels are intended, and how was their tissue acquired?
4. Which methodological warnings must be resolved before quantitative comparison?

The reference panels currently appear as `declared_not_materialized`. This is intentional. The
configuration is ready for the GTEx and Tabula Sapiens adapters, but the current release has not
downloaded or harmonized those matrices.

## What this prevents

The workflow prevents several quiet analytical errors:

- Calling adjacent tissue healthy.
- Treating absent comorbidity metadata as evidence that no comorbidity existed.
- Counting multiple samples from one person as independent donors.
- Presenting independently processed disease and reference matrices as an unconfounded contrast.
- Losing the exact query and retrieval time used to construct a cohort.

## Next extension

The GTEx adapter will materialize donor/sample metadata before expression data. That enables age and
sex overlap plots, source-specific missingness, and eligibility-context review before any molecular
comparison. The Tabula Sapiens adapter will then add donor-aware cell-type summaries while keeping
organ-donor context and study effects visible.

