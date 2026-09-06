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

`comparability` records intended stratification, minimum donor counts, handling of missing metadata,
and whether source effects must remain visible.
