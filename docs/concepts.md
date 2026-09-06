# Cohorts and reference context

## A reference population answers a particular question

A reference panel is not automatically an estimate of an ideal healthy population. Its meaning
depends on donor eligibility, age and comorbidities, tissue procurement, processing, assay, and
metadata completeness. Open Cohort Factory therefore models reference context separately from
tissue identity.

| Context | Interpretation |
| --- | --- |
| `adjacent_non_tumor` | Tissue collected from a person with cancer, near or separate from a tumor |
| `postmortem_reference` | Tissue obtained after death under source-specific eligibility criteria |
| `organ_donor` | Tissue collected through an organ-donation process |
| `surgical_non_diseased` | Tissue obtained during surgery and considered non-diseased for the study |
| `healthy_volunteer` | Tissue or cells collected from a specifically enrolled healthy volunteer |
| `disease_control` | Non-index-disease comparator with another known condition |
| `unknown_or_incompletely_characterized` | Acquisition or health context cannot be established |

These labels describe provenance. They do not guarantee molecular comparability.

## Missing does not mean absent

If a public dataset lacks comorbidity, medication, smoking, or other exposure metadata, the software
reports the field as unknown. It must not infer that an unrecorded condition or exposure was absent.

## Donor-aware analysis

Multiple samples may come from one donor, and single-cell datasets contain many cells per donor.
The manifest retains those relationships. Descriptive reports show both sample and donor counts;
future single-cell analyses will use donors, not cells, as the independent inferential units.

## Cross-study comparisons

Source, assay, tissue processing, ischemic time, and population composition can all be confounded
with disease status. The project will not present a naïve TCGA-versus-GTEx comparison as biological
differential expression. Reference sources remain separately visible unless a documented,
assay-aware method justifies integration.

