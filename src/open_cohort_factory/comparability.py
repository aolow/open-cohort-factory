"""Donor-level demographic and specimen-context comparisons."""

from __future__ import annotations

from collections import Counter
from typing import Any

from .models import SampleRecord


def _unique_donors(samples: list[SampleRecord]) -> list[SampleRecord]:
    donors: dict[tuple[str, str, str], SampleRecord] = {}
    for sample in samples:
        key = (sample.source.value, sample.cohort_name, sample.case_id)
        donors.setdefault(key, sample)
    return list(donors.values())


def _proportion(counts: Counter[str], category: str) -> float | None:
    total = sum(counts.values())
    return round(counts[category] / total, 3) if total else None


def compare_populations(samples: list[SampleRecord]) -> list[dict[str, Any]]:
    """Compare each materialized reference panel with the disease cohort."""
    donors = _unique_donors(samples)
    disease = [sample for sample in donors if sample.cohort_role == "disease"]
    disease_age = Counter(sample.age_bracket or "unknown" for sample in disease)
    disease_sex = Counter(sample.sex_at_birth or "unknown" for sample in disease)
    comparisons: list[dict[str, Any]] = []

    reference_names = sorted(
        {sample.cohort_name for sample in donors if sample.cohort_role == "reference"}
    )
    for name in reference_names:
        reference = [sample for sample in donors if sample.cohort_name == name]
        reference_age = Counter(sample.age_bracket or "unknown" for sample in reference)
        reference_sex = Counter(sample.sex_at_birth or "unknown" for sample in reference)
        disease_brackets = {key for key in disease_age if key != "unknown"}
        reference_brackets = {key for key in reference_age if key != "unknown"}
        uncovered = sorted(disease_brackets - reference_brackets)
        disease_female = _proportion(disease_sex, "female")
        reference_female = _proportion(reference_sex, "female")
        sex_difference = (
            round(abs(disease_female - reference_female), 3)
            if disease_female is not None and reference_female is not None
            else None
        )
        pathology = Counter(
            category for sample in reference for category in sample.pathology_categories_present
        )

        if uncovered:
            decision = "context_only_for_full_cohort"
        elif sex_difference is not None and sex_difference >= 0.1:
            decision = "requires_stratification_or_weighting"
        else:
            decision = "demographically_comparable_with_caveats"

        comparisons.append(
            {
                "reference_panel": name,
                "disease_donors": len(disease),
                "reference_donors": len(reference),
                "disease_age_brackets": dict(disease_age),
                "reference_age_brackets": dict(reference_age),
                "uncovered_disease_age_brackets": uncovered,
                "female_proportion": {
                    "disease": disease_female,
                    "reference": reference_female,
                    "absolute_difference": sex_difference,
                },
                "top_reference_pathology_categories": [
                    {"category": category, "samples": count}
                    for category, count in pathology.most_common(10)
                ],
                "decision": decision,
                "recommended_actions": _recommendations(uncovered, sex_difference, bool(pathology)),
            }
        )
    return comparisons


def _recommendations(
    uncovered_age_brackets: list[str], sex_difference: float | None, has_pathology: bool
) -> list[str]:
    actions = ["Keep source-specific results separate in all descriptive summaries."]
    if uncovered_age_brackets:
        actions.append(
            "Restrict the disease cohort to reference-covered ages or report uncovered ages "
            "as a separate sensitivity stratum."
        )
    if sex_difference is not None and sex_difference >= 0.1:
        actions.append("Stratify or weight by sex before comparing population summaries.")
    if has_pathology:
        actions.append(
            "Run sensitivity analyses using GTEx specimen pathology annotations; do not equate "
            "unrecorded pathology with its absence."
        )
    actions.append(
        "Do not run cross-study differential expression without assay-aware harmonization."
    )
    return actions
