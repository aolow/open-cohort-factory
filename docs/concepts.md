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

Missingness can itself have structure. Race, ethnicity, procurement timing, RNA integrity, death
classification, and pathology may be available in one source but absent in another. In that case,
the covariate is confounded with source at the design level. Restricting to complete cases or adding
the field to a regression model cannot recover information that was never collected in one arm.
The report also labels covariate roles. Post-baseline outcomes such as vital status should not be
treated interchangeably with baseline adjustment variables, because conditioning on an outcome of
disease or treatment can introduce rather than remove bias.

## Donor-aware analysis

Multiple samples may come from one donor, and single-cell datasets contain many cells per donor.
The manifest retains those relationships. Descriptive reports show both sample and donor counts;
future single-cell analyses will use donors, not cells, as the independent inferential units.

## Uniform processing is not population comparability

Source, assay, tissue processing, ischemic time, and population composition can all be confounded
with disease status. The UCSC Xena Toil matrix reduces processing differences by applying a common
RNA-seq pipeline to TCGA, TARGET, and GTEx. That makes descriptive expression contrasts useful, but
does not make the contributing populations exchangeable. Open Cohort Factory therefore selects and
audits the populations before requesting their expression values and keeps reference sources visible.

## Matching is a design aid, not harmonization

The initial matching method creates deterministic 1:1 donor pairs within exact public age-bracket
and sex strata. Selection within an oversubscribed stratum uses a recorded seed and stable hash, so
the result is reproducible without implying that identifier order is biologically meaningful.

Matching improves balance only for included variables. It does not resolve tissue procurement,
postmortem effects, assay processing, ancestry, smoking, treatment, comorbidity, or unmeasured
differences. Unavailable metadata cannot be balanced. Consequently, a matched donor export is an
auditable analysis population—not permission to perform naïve cross-study differential expression.

When tumor and adjacent non-tumor tissue come from the same participant, the participant identifier
takes precedence over demographic matching. This controls participant-level differences but does
not make adjacent tissue healthy: field effects, occult alteration, and cancer-associated systemic
effects remain plausible.

## Uncertainty follows the donor design

Expression summaries use donor medians so multiple biospecimens do not become independent people.
The report calculates deterministic 95% bootstrap intervals by resampling donors for unmatched
contrasts and resampling complete donor pairs together for demographic matches and within-participant
tumor/adjacent comparisons. These intervals describe sampling uncertainty in the assembled public
cohorts; they do not account for every unmeasured biological or technical source effect.
