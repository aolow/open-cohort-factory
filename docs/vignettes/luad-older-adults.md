# Older adults with lung adenocarcinoma

This vignette builds a public TCGA lung adenocarcinoma cohort restricted to adults aged 60 or older,
materializes bulk postmortem lung metadata from GTEx, adds cell-type-level organ-donor context from
Tabula Sapiens, and retrieves a small expression slice from the UCSC Xena Toil recompute.

The purpose is to show that technical harmonization and population comparability are separate
requirements: Xena supplies consistently processed values, while cohort construction determines
which donors belong in each descriptive comparison.

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
  - name: TCGA-LUAD adjacent non-tumor lung
    source: gdc
    tissue: Bronchus and lung
    context: adjacent_non_tumor
    sample_types: [Solid Tissue Normal]

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

expression:
  source: xena_toil
  genes: [EPCAM, CEACAM5, MSLN]
  dataset: TcgaTargetGtex_rsem_gene_tpm
  transform: log2_tpm_plus_1
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
├── expression.tsv
├── matched_donors.tsv
├── matching.json
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

All three reference panels appear as `materialized`. Xena expression is attached only to eligible
TCGA and GTEx identifiers; Tabula Sapiens remains cellular context and is not treated as though it
shares the bulk-expression matrix.

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

## Build a matched analysis population

The example enables deterministic 1:1 exact matching on public age bracket and sex:

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

In the September 6, 2026 verification run, this yielded 171 donor pairs:

| Age bracket | Sex | Matched pairs |
| --- | --- | ---: |
| 60–69 | Female | 66 |
| 60–69 | Male | 83 |
| 70–79 | Female | 4 |
| 70–79 | Male | 18 |

The 342 selected donor rows are written to `matched_donors.tsv`. The complete audit in
`matching.json` retains 242 donors that were not selected and records whether each was excluded due
to a missing opposite-population stratum, surplus within a stratum, or missing matching metadata.

This matched population is suitable for controlled descriptive work. Exact demographic balance does
not remove postmortem, procurement, residual study, or unmeasured clinical differences.

## Compare eligible and expression-complete populations

In the verified Xena-backed run, 738 unique cohort sample identifiers were requested and 474 were
present in the Toil matrix: 357 TCGA-LUAD tumors, 42 TCGA adjacent samples, and 75 GTEx lung
samples. This incomplete overlap is retained in provenance rather than silently ignored; the current
source releases and the older Toil compendium do not contain identical sample sets.

Of the 171 GTEx demographic donor pairs, 60 had Xena measurements for both members. Matched
expression summaries therefore use those 60 complete pairs, with donor medians as the analytical
unit. Separately, 42 TCGA participants had both primary-tumor and adjacent-tissue expression:

| Gene | All GTEx | Matched GTEx | All adjacent | Within-participant adjacent |
| --- | ---: | ---: | ---: | ---: |
| CEACAM5 | 4.8284 | 4.6023 | 3.4664 | 2.1239 |
| EPCAM | 3.3376 | 3.0989 | 1.4166 | 1.2072 |
| MSLN | 0.9434 | 1.1113 | -0.3796 | -0.6131 |

Values are differences of medians on the `log2(TPM + 1)` scale. They are descriptive effect
summaries, not claims that postmortem GTEx lung represents a universally healthy counterfactual.
The report adds donor-level IQRs and 95% bootstrap intervals; paired populations are resampled as
pairs. In the verified run, the interval for MSLN crossed zero in the matched GTEx and both adjacent
analyses, while CEACAM5 and EPCAM remained positive. This distinguishes an unstable reference-sensitive
signal from a consistently directed one.

## Audit hidden covariates

The covariate landscape shows both distributions and whether each variable was collected. In this
example, age and sex are broadly available, but race, ethnicity, and vital status are available in
GDC and absent from the selected GTEx and Tabula Sapiens metadata. Conversely, ischemic time, RIN,
Hardy death classification, and detailed specimen pathology are available for GTEx but absent from
the GDC cohort representation. These are not ordinary random missing values: availability is tied to
source and therefore to disease/reference status.

The practical conclusion is that matching age and sex improves the comparison without making it
fully adjusted. Procurement and specimen-quality effects remain entangled with the definition of
the reference population and must stay visible in interpretation.

## Choose a reference for a question, not in the abstract

The reference-fitness panel turns the audit into a decision aid without naming one universal winner.
In this example, adjacent TCGA lung is the natural choice for a within-participant sensitivity
analysis; matched GTEx is useful for adult tissue-expression context with procurement caveats; and
Tabula Sapiens is useful for cell-type attribution but not as a bulk-expression or population-level
counterfactual. None supports an unqualified claim about expression in a universally healthy person.

## Screen expression beyond the tissue of origin

The pan-tissue option queried 7,425 GTEx samples across 30 primary-tissue groups in the verified
Toil run. At a threshold of 1 TPM, at least half of donors were positive in 7 tissues for CEACAM5,
21 tissues for EPCAM, and 7 tissues for MSLN.

The highest normal-tissue signals materially change interpretation. CEACAM5 was highest in vagina,
colon, small intestine, and salivary gland; EPCAM was broadly expressed, led by small intestine,
thyroid, pituitary, and salivary gland; MSLN was highest in lung, followed by fallopian tube and
salivary gland. The five-donor fallopian-tube estimate is displayed with its donor count rather than
given the same evidentiary weight as tissues represented by hundreds of donors.

This screen is an expression-context flag, not a toxicity prediction. Bulk RNA does not establish
protein abundance, accessibility, essentiality, or the cell type carrying the signal; those questions
motivate the next cell-type-attribution phase.

With `cell_type_attribution: true`, the build retrieves the same genes from the declared Tabula
Sapiens lung panel. It reports detection and raw-count summaries at donor × cell-type resolution,
then ranks cell types using donor medians. This distinguishes a signal distributed through epithelial
populations from one concentrated in a rarer immune or stromal compartment.

The attribution is nested in a collapsed report section. It answers “which annotated cells may
contribute to this tissue signal?”—not whether a protein is on the cell surface, accessible to a
therapeutic modality, or causally responsible for toxicity. Cell labels, dissociation, sampling,
dropout, and limited donor counts remain important limitations.

In the verified lung run, EPCAM was most consistently detected in alveolar type 2, multiciliated,
club, and goblet epithelial cells. CEACAM5 was led by goblet and club cells. MSLN was led by
mesothelial cells, but that estimate came from only 11 cells across two donors; the report retains
both counts so an apparently strong signal cannot hide its limited support.

## Add cellular context without inflating the evidence

The same build queries the tissue-specific Tabula Sapiens lung dataset through CELLxGENE Census.
For the configured age of 60 or older, the verified run returned:

- 36,364 cells
- 89 donor × cell-type aggregate records
- Two eligible donors: one female and one male
- Cell ontology identifiers retained alongside labels

The most abundant captured populations included alveolar type 2 cells, macrophages, capillary
endothelial cells, basal cells, monocytes, club cells, and T cells. The report shows donor coverage
beside every cell count; a population with thousands of cells from two people is still evidence from
two people.

The comparability decision is therefore `insufficient_reference_donors`. Tabula Sapiens can inform
which cell types may contribute to a bulk lung signal, but this subset cannot support population-level
inference about older adults with lung cancer. Even though the matcher can construct two exact
demographic pairs, those pairs are retained only as an auditable demonstration—not as an adequately
powered inferential cohort.

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

The sensitivity result is the point of the workflow: the apparent tumor-reference difference changes
with reference definition. The report does not select whichever comparator produces the most
favorable effect; it presents the acquisition context, overlap, and attrition needed to interpret
each estimate.
