"""Deterministic, auditable donor-level cohort matching."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from hashlib import sha256
from typing import Any

from .models import MatchingSpec, SampleRecord


def _unique_donors(samples: list[SampleRecord]) -> list[SampleRecord]:
    donors: dict[tuple[str, str, str], SampleRecord] = {}
    for sample in samples:
        key = (sample.source.value, sample.cohort_name, sample.case_id)
        donors.setdefault(key, sample)
    return list(donors.values())


def _stratum(sample: SampleRecord, fields: Sequence[str]) -> tuple[str, ...] | None:
    values = [getattr(sample, field) for field in fields]
    if any(value is None or value == "unknown" for value in values):
        return None
    return tuple(str(value) for value in values)


def _rank(sample: SampleRecord, seed: int) -> str:
    value = f"{seed}|{sample.source.value}|{sample.cohort_name}|{sample.case_id}"
    return sha256(value.encode()).hexdigest()


def exact_match(samples: list[SampleRecord], spec: MatchingSpec) -> list[dict[str, Any]]:
    """Create independent 1:1 exact matches for each materialized reference panel."""
    if not spec.enabled:
        return []
    donors = _unique_donors(samples)
    disease = [sample for sample in donors if sample.cohort_role == "disease"]
    reference_names = sorted(
        {sample.cohort_name for sample in donors if sample.cohort_role == "reference"}
    )
    results: list[dict[str, Any]] = []
    for reference_name in reference_names:
        reference = [sample for sample in donors if sample.cohort_name == reference_name]
        disease_groups = _group_by_stratum(disease, spec.variables)
        reference_groups = _group_by_stratum(reference, spec.variables)
        selected: list[dict[str, Any]] = []
        excluded: list[dict[str, Any]] = []
        pair_number = 0

        for stratum in sorted(set(disease_groups) | set(reference_groups)):
            disease_group = sorted(disease_groups[stratum], key=lambda row: _rank(row, spec.seed))
            reference_group = sorted(
                reference_groups[stratum], key=lambda row: _rank(row, spec.seed)
            )
            matched_count = min(len(disease_group), len(reference_group))
            for disease_donor, reference_donor in zip(
                disease_group[:matched_count], reference_group[:matched_count], strict=True
            ):
                pair_number += 1
                pair_id = f"{reference_name}-pair-{pair_number:04d}"
                selected.extend(
                    [
                        _selection_row(pair_id, disease_donor, stratum, spec.variables),
                        _selection_row(pair_id, reference_donor, stratum, spec.variables),
                    ]
                )
            excluded.extend(
                _exclusion_rows(disease_group[matched_count:], stratum, "disease", reference_group)
            )
            excluded.extend(
                _exclusion_rows(
                    reference_group[matched_count:], stratum, "reference", disease_group
                )
            )

        for sample in disease + reference:
            if _stratum(sample, spec.variables) is None:
                excluded.append(
                    {
                        "source": sample.source.value,
                        "cohort_role": sample.cohort_role,
                        "cohort_name": sample.cohort_name,
                        "donor_id": sample.case_id,
                        "reason": "missing_matching_field",
                    }
                )

        results.append(
            {
                "reference_panel": reference_name,
                "method": spec.method,
                "variables": spec.variables,
                "ratio": spec.ratio,
                "seed": spec.seed,
                "matched_pairs": pair_number,
                "selected_donors": selected,
                "excluded_donors": excluded,
                "matched_strata": _matched_strata(selected, spec.variables),
            }
        )
    return results


def _group_by_stratum(
    samples: list[SampleRecord], fields: Sequence[str]
) -> defaultdict[tuple[str, ...], list[SampleRecord]]:
    groups: defaultdict[tuple[str, ...], list[SampleRecord]] = defaultdict(list)
    for sample in samples:
        stratum = _stratum(sample, fields)
        if stratum is not None:
            groups[stratum].append(sample)
    return groups


def _selection_row(
    pair_id: str, sample: SampleRecord, stratum: tuple[str, ...], fields: Sequence[str]
) -> dict[str, Any]:
    return {
        "pair_id": pair_id,
        "source": sample.source.value,
        "cohort_role": sample.cohort_role,
        "cohort_name": sample.cohort_name,
        "donor_id": sample.case_id,
        "stratum": dict(zip(fields, stratum, strict=True)),
    }


def _exclusion_rows(
    samples: list[SampleRecord],
    stratum: tuple[str, ...],
    role: str,
    opposite_group: list[SampleRecord],
) -> list[dict[str, Any]]:
    reason = "no_opposite_donors_in_stratum" if not opposite_group else "surplus_in_stratum"
    return [
        {
            "source": sample.source.value,
            "cohort_role": role,
            "cohort_name": sample.cohort_name,
            "donor_id": sample.case_id,
            "stratum": list(stratum),
            "reason": reason,
        }
        for sample in samples
    ]


def _matched_strata(selected: list[dict[str, Any]], fields: Sequence[str]) -> list[dict[str, Any]]:
    counts: dict[tuple[str, ...], int] = defaultdict(int)
    for row in selected:
        if row["cohort_role"] == "disease":
            key = tuple(row["stratum"][field] for field in fields)
            counts[key] += 1
    return [
        {"stratum": dict(zip(fields, key, strict=True)), "matched_pairs": count}
        for key, count in sorted(counts.items())
    ]
