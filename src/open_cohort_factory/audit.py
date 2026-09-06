"""Cohort metadata and comparability audits."""

from __future__ import annotations

import statistics
from collections import Counter
from typing import Any

from .comparability import compare_populations
from .models import ProjectSpec, SampleRecord


def summarize(samples: list[SampleRecord], spec: ProjectSpec) -> dict[str, Any]:
    disease = [sample for sample in samples if sample.cohort_role == "disease"]
    ages = [sample.age_at_index for sample in disease if sample.age_at_index is not None]
    missing = Counter(field for sample in samples for field in sample.metadata_missing)
    population_rows = []
    populations = sorted({(sample.cohort_role, sample.cohort_name) for sample in samples})
    for role, name in populations:
        members = [
            sample for sample in samples if (sample.cohort_role, sample.cohort_name) == (role, name)
        ]
        donor_members = list({sample.case_id: sample for sample in members}.values())
        population_rows.append(
            {
                "role": role,
                "name": name,
                "source": members[0].source.value,
                "samples": len(members),
                "donors": len(donor_members),
                "cells": sum(sample.cell_count or 0 for sample in members) or None,
                "age_brackets": dict(
                    Counter(sample.age_bracket or "unknown" for sample in donor_members)
                ),
                "sex": dict(Counter(sample.sex_at_birth or "unknown" for sample in donor_members)),
            }
        )
    return {
        "sample_count": len(samples),
        "donor_count": len({(sample.source, sample.case_id) for sample in samples}),
        "project_counts": dict(Counter(sample.project_id for sample in samples)),
        "sample_type_counts": dict(Counter(sample.sample_type or "unknown" for sample in samples)),
        "sex_counts": dict(Counter(sample.sex_at_birth or "unknown" for sample in samples)),
        "age": {
            "known": len(ages),
            "minimum": min(ages) if ages else None,
            "maximum": max(ages) if ages else None,
            "mean": round(sum(ages) / len(ages), 1) if ages else None,
        },
        "missing_fields": dict(missing),
        "populations": population_rows,
        "comparability": compare_populations(
            samples, minimum_reference_donors=spec.comparability.minimum_donors_per_group
        ),
        "cell_type_summaries": _cell_type_summaries(samples),
        "covariate_landscape": _covariate_landscape(samples),
        "reference_panels": [
            {
                "name": panel.name,
                "source": panel.source.value,
                "context": panel.context.value,
                "tissue": panel.tissue,
                "resolution": panel.resolution,
                "status": (
                    "materialized"
                    if any(sample.cohort_name == panel.name for sample in samples)
                    else "declared_not_materialized"
                ),
            }
            for panel in spec.reference_panels
        ],
    }


COVARIATES = [
    ("age_bracket", "Age bracket", "categorical", "demographic"),
    ("sex_at_birth", "Sex at birth", "categorical", "demographic"),
    ("race", "Race", "categorical", "social-demographic proxy"),
    ("ethnicity", "Ethnicity", "categorical", "social-demographic proxy"),
    ("vital_status", "Vital status", "categorical", "post-baseline outcome"),
    ("assay_type", "Assay type", "categorical", "technical"),
    ("ischemic_time_minutes", "Ischemic time", "numeric", "procurement"),
    ("rin", "RNA integrity (RIN)", "numeric", "specimen quality"),
    ("hardy_scale", "Death classification", "categorical", "procurement"),
    ("pathology_notes", "Specimen pathology", "categorical", "specimen condition"),
]


def _known_covariate(value: Any) -> bool:
    if value is None or value == "":
        return False
    if isinstance(value, str):
        return value.strip().casefold() not in {
            "unknown",
            "not reported",
            "not available",
            "not specified",
        }
    return True


def _covariate_landscape(samples: list[SampleRecord]) -> dict[str, Any]:
    populations = sorted({(sample.cohort_role, sample.cohort_name) for sample in samples})
    donor_groups: dict[str, list[SampleRecord]] = {}
    for role, name in populations:
        members = [
            sample for sample in samples if (sample.cohort_role, sample.cohort_name) == (role, name)
        ]
        donor_groups[name] = list({sample.case_id: sample for sample in members}.values())

    fields: list[dict[str, Any]] = []
    for field, label, kind, covariate_role in COVARIATES:
        coverage = []
        known_rates = []
        category_distributions: list[dict[str, float]] = []
        for population, donors in donor_groups.items():
            values = [getattr(donor, field) for donor in donors]
            known = [value for value in values if _known_covariate(value)]
            rate = len(known) / len(donors) if donors else 0.0
            known_rates.append(rate)
            detail: dict[str, Any] = {}
            if known and kind == "numeric":
                numeric = [float(value) for value in known]
                detail = {
                    "median": round(statistics.median(numeric), 2),
                    "minimum": round(min(numeric), 2),
                    "maximum": round(max(numeric), 2),
                }
            elif known:
                counts = Counter(str(value) for value in known)
                if len(known) >= 5:
                    category_distributions.append(
                        {category: count / len(known) for category, count in counts.items()}
                    )
                detail = {
                    "top_values": [
                        {"value": str(value), "donors": count}
                        for value, count in counts.most_common(3)
                    ]
                }
            coverage.append(
                {
                    "population": population,
                    "known": len(known),
                    "total": len(donors),
                    "percent": round(rate * 100),
                    "detail": detail,
                }
            )
        spread = max(known_rates, default=0) - min(known_rates, default=0)
        distribution_distance = _maximum_distribution_distance(category_distributions)
        if max(known_rates, default=0) == 0:
            availability_status = "unavailable"
        elif spread >= 0.5 and max(known_rates, default=0) >= 0.5:
            availability_status = "source_confounded"
        else:
            availability_status = "partially_comparable"
        fields.append(
            {
                "field": field,
                "label": label,
                "kind": kind,
                "role": covariate_role,
                "coverage": coverage,
                "source_confounded_availability": (
                    spread >= 0.5 and max(known_rates, default=0) >= 0.5
                ),
                "availability_status": availability_status,
                "distribution_shift": distribution_distance >= 0.2,
                "maximum_distribution_distance": round(distribution_distance, 3),
            }
        )
    return {
        "populations": list(donor_groups),
        "fields": fields,
        "interpretation": (
            "A covariate measured primarily in one source is structurally confounded with source; "
            "missingness cannot be repaired by ordinary adjustment."
        ),
    }


def _maximum_distribution_distance(distributions: list[dict[str, float]]) -> float:
    maximum = 0.0
    for index, left in enumerate(distributions):
        for right in distributions[index + 1 :]:
            categories = left.keys() | right.keys()
            distance = 0.5 * sum(abs(left.get(key, 0) - right.get(key, 0)) for key in categories)
            maximum = max(maximum, distance)
    return maximum


def methodological_warnings(spec: ProjectSpec, samples: list[SampleRecord]) -> list[str]:
    warnings = [
        "Reference samples are context-specific; the report does not treat any panel as a "
        "universal normal baseline.",
        "Cross-study expression comparisons require assay-aware harmonization and must not be "
        "interpreted as unconfounded differential expression.",
        "Missing clinical metadata means unknown, not absence of a condition or exposure.",
    ]
    if any(
        sample.reference_context is not None
        and sample.reference_context.value == "adjacent_non_tumor"
        for sample in samples
    ):
        warnings.append(
            "Adjacent non-tumor samples come from people with cancer and are not equivalent to "
            "healthy tissue."
        )
    if len({panel.source for panel in spec.reference_panels}) > 1:
        warnings.append("Reference-source effects must remain visible in all combined summaries.")
    if any(sample.source.value == "gtex" for sample in samples):
        warnings.append(
            "Public GTEx age is grouped into decade brackets; demographic overlap is therefore "
            "coarser than the exact-age disease-cohort summary."
        )
        warnings.append(
            "GTEx pathology and death-related metadata describe the donated specimen but do not "
            "establish the absence of every comorbidity or exposure."
        )
    return warnings


def _cell_type_summaries(samples: list[SampleRecord]) -> list[dict[str, Any]]:
    panels = sorted({sample.cohort_name for sample in samples if sample.cell_type is not None})
    results: list[dict[str, Any]] = []
    for panel in panels:
        members = [sample for sample in samples if sample.cohort_name == panel]
        cell_types = sorted(
            {sample.cell_type for sample in members if sample.cell_type is not None}
        )
        rows: list[dict[str, Any]] = []
        for cell_type in cell_types:
            typed = [sample for sample in members if sample.cell_type == cell_type]
            rows.append(
                {
                    "cell_type": cell_type,
                    "ontology_term_id": typed[0].cell_type_ontology_term_id,
                    "cells": sum(sample.cell_count or 0 for sample in typed),
                    "donors": len({sample.case_id for sample in typed}),
                }
            )
        results.append(
            {
                "reference_panel": panel,
                "total_cells": sum(row["cells"] for row in rows),
                "total_donors": len({sample.case_id for sample in members}),
                "cell_types": sorted(rows, key=lambda row: row["cells"], reverse=True),
            }
        )
    return results
