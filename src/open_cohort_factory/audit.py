"""Cohort metadata and comparability audits."""

from __future__ import annotations

from collections import Counter
from typing import Any

from .models import ProjectSpec, SampleRecord


def summarize(samples: list[SampleRecord], spec: ProjectSpec) -> dict[str, Any]:
    ages = [sample.age_at_index for sample in samples if sample.age_at_index is not None]
    missing = Counter(field for sample in samples for field in sample.metadata_missing)
    return {
        "sample_count": len(samples),
        "donor_count": len({sample.case_id for sample in samples}),
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
        "reference_panels": [
            {
                "name": panel.name,
                "source": panel.source.value,
                "context": panel.context.value,
                "tissue": panel.tissue,
                "resolution": panel.resolution,
                "status": "declared_not_materialized",
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
    if any(sample.reference_context for sample in samples):
        warnings.append(
            "Adjacent non-tumor samples come from people with cancer and are not equivalent to "
            "healthy tissue."
        )
    if len({panel.source for panel in spec.reference_panels}) > 1:
        warnings.append("Reference-source effects must remain visible in all combined summaries.")
    return warnings
