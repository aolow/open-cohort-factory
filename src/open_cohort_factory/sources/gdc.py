"""NCI Genomic Data Commons adapter."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from typing import Any

import httpx

from ..models import (
    DataSource,
    DiseaseCohortSpec,
    ProvenanceRecord,
    ReferenceContext,
    SampleRecord,
)

GDC_CASES_ENDPOINT = "https://api.gdc.cancer.gov/cases"
GDC_FIELDS = [
    "case_id",
    "project.project_id",
    "primary_site",
    "demographic.age_at_index",
    "demographic.sex_at_birth",
    "demographic.gender",
    "demographic.race",
    "demographic.ethnicity",
    "demographic.vital_status",
    "diagnoses.age_at_diagnosis",
    "diagnoses.tissue_or_organ_of_origin",
    "samples.sample_id",
    "samples.sample_type",
]


def _clause(field: str, values: list[Any], operator: str = "in") -> dict[str, Any]:
    return {"op": operator, "content": {"field": field, "value": values}}


def build_filters(spec: DiseaseCohortSpec) -> dict[str, Any]:
    clauses: list[dict[str, Any]] = [_clause("project.project_id", spec.projects)]
    if spec.primary_sites:
        clauses.append(_clause("primary_site", spec.primary_sites))
    if spec.sample_types:
        clauses.append(_clause("samples.sample_type", spec.sample_types))
    if spec.sex_at_birth:
        clauses.append(_clause("demographic.sex_at_birth", spec.sex_at_birth))
    if spec.age.minimum is not None:
        clauses.append(_clause("demographic.age_at_index", [spec.age.minimum], ">="))
    if spec.age.maximum is not None:
        clauses.append(_clause("demographic.age_at_index", [spec.age.maximum], "<="))
    return {"op": "and", "content": clauses}


def _first(items: list[dict[str, Any]] | None, key: str) -> Any:
    return items[0].get(key) if items else None


def _reference_context(sample_type: str | None) -> ReferenceContext | None:
    if sample_type and "normal" in sample_type.lower():
        return ReferenceContext.ADJACENT_NON_TUMOR
    return None


class GDCClient:
    """Retrieve public case and sample metadata from GDC."""

    def __init__(self, client: httpx.Client | None = None) -> None:
        self._client = client or httpx.Client(timeout=60, follow_redirects=True)

    def fetch(
        self, spec: DiseaseCohortSpec, page_size: int = 500
    ) -> tuple[list[SampleRecord], ProvenanceRecord]:
        filters = build_filters(spec)
        params: dict[str, Any] = {
            "filters": json.dumps(filters),
            "format": "JSON",
            "fields": ",".join(GDC_FIELDS),
            "expand": "project,demographic,diagnoses,samples",
            "size": page_size,
        }
        hits: list[dict[str, Any]] = []
        offset = 0
        while True:
            params["from"] = offset
            response = self._client.get(GDC_CASES_ENDPOINT, params=params)
            response.raise_for_status()
            data = response.json().get("data", {})
            page = data.get("hits", [])
            hits.extend(page)
            total = data.get("pagination", {}).get("total", len(hits))
            offset += len(page)
            if not page or offset >= total:
                break
        samples = [sample for hit in hits for sample in self._normalize_case(hit, spec)]
        provenance = ProvenanceRecord(
            source=DataSource.GDC,
            endpoint=GDC_CASES_ENDPOINT,
            retrieved_at=datetime.now(UTC).isoformat(),
            query=filters,
            citation_url="https://gdc.cancer.gov/developers/gdc-application-programming-interface-api",
        )
        return samples, provenance

    @staticmethod
    def _normalize_case(hit: dict[str, Any], spec: DiseaseCohortSpec) -> list[SampleRecord]:
        demographic = hit.get("demographic") or {}
        diagnoses = hit.get("diagnoses") or []
        project_id = (hit.get("project") or {}).get("project_id", "unknown")
        records: list[SampleRecord] = []
        for sample in hit.get("samples") or []:
            sample_type = sample.get("sample_type")
            if spec.sample_types and sample_type not in spec.sample_types:
                continue
            record = SampleRecord(
                source=DataSource.GDC,
                cohort_role="disease",
                cohort_name="disease_cohort",
                project_id=project_id,
                case_id=hit["case_id"],
                sample_id=sample["sample_id"],
                sample_type=sample_type,
                primary_site=hit.get("primary_site"),
                tissue_or_organ_of_origin=_first(diagnoses, "tissue_or_organ_of_origin"),
                age_at_index=demographic.get("age_at_index"),
                age_bracket=_age_bracket(demographic.get("age_at_index")),
                age_at_diagnosis_days=_first(diagnoses, "age_at_diagnosis"),
                sex_at_birth=demographic.get("sex_at_birth") or demographic.get("gender"),
                race=demographic.get("race"),
                ethnicity=demographic.get("ethnicity"),
                vital_status=demographic.get("vital_status"),
                reference_context=_reference_context(sample_type),
                acquisition_context="cancer study participant",
                source_payload=hit,
            )
            required = ["age_at_index", "sex_at_birth", "primary_site", "sample_type"]
            record.metadata_missing = [name for name in required if getattr(record, name) is None]
            records.append(record)
        return records


def _age_bracket(age: int | None) -> str | None:
    if age is None:
        return None
    lower = age // 10 * 10
    return f"{lower}-{lower + 9}"
