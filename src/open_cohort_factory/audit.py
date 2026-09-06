"""Cohort metadata and comparability audits."""

from __future__ import annotations

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
