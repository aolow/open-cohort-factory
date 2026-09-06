# Command reference

## Validate

```bash
uv run cohort-factory validate PATH_TO_SPEC.yaml
```

Parses the YAML specification and validates all fields without accessing the network or writing
outputs.

## Build

```bash
uv run cohort-factory build PATH_TO_SPEC.yaml --output OUTPUT_DIRECTORY
```

Retrieves matching public GDC disease-cohort and GTEx reference-panel metadata, performs donor-level
comparability checks, and writes a manifest, summary, and HTML report. Existing files with the same
names in the selected output directory are replaced.

## Specification sections

`disease_cohort` defines the currently materialized GDC population. Supported filters include GDC
project identifiers, primary site, age range, sex at birth, and sample type.

`reference_panels` declares one or more contextual comparator populations. Each panel requires a
source, tissue, acquisition context, resolution, and optional age limits and notes. GTEx panels with
bulk resolution are currently materialized through the GTEx Portal V2 API using RNA-seq samples.
Tabula Sapiens panels with cell-type resolution are queried through the pinned CELLxGENE Census LTS
release and emitted as donor × cell-type aggregate records. GDC panels may declare `projects` and
must declare `sample_types`; this supports adjacent non-tumor tissue without relabeling it healthy.

`comparability` records intended stratification, minimum donor counts, handling of missing metadata,
and whether source effects must remain visible.

`matching` controls optional deterministic donor matching. The current implementation supports
1:1 exact matching on `age_bracket`, `sex_at_birth`, or both. `seed` controls reproducible selection
when a stratum contains more eligible donors than needed. Adjacent non-tumor panels are linked to
tumor samples from the same GDC participant instead of being randomly paired within demographic
strata.

```yaml
matching:
  enabled: true
  method: exact
  variables:
    - age_bracket
    - sex_at_birth
  ratio: 1
  seed: 2026
```

Selected pairs are written to `matched_donors.tsv`. `matching.json` also retains excluded donors and
distinguishes missing matching fields, absent opposite-population strata, and surplus donors.

`expression` optionally requests a small gene slice from the uniformly processed UCSC Xena Toil
matrix after cohort construction and matching. TCGA, TARGET, and GTEx labels are never collapsed;
only GDC/TCGA and GTEx samples explicitly present in the constructed cohort are queried.

```yaml
expression:
  source: xena_toil
  genes: [EPCAM, CEACAM5, MSLN]
  dataset: TcgaTargetGtex_rsem_gene_tpm
  transform: log2_tpm_plus_1
```

Xena values arrive as `log2(TPM + 0.001)` and are re-expressed as `log2(TPM + 1)`. The manifest
records the dataset, genes, requested and returned sample counts, source-specific coverage, input
scale, output transform, retrieval time, reference build, and annotation version. Matched summaries
retain only pairs for which both donors have Xena measurements and use donor medians when a donor
has multiple samples.
