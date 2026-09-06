"""GTEx Portal V2 API adapter for reference-tissue metadata."""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Any

import httpx

from ..models import DataSource, ProvenanceRecord, ReferencePanelSpec, SampleRecord

GTEX_SAMPLE_ENDPOINT = "https://gtexportal.org/api/v2/dataset/sample"
GTEX_DATASET_ID = "gtex_v10"
GTEX_CITATION_URL = "https://gtexportal.org/home/downloads/adult-gtex/overview"


def age_brackets(minimum: int | None, maximum: int | None) -> list[str]:
    """Return public GTEx decade brackets that overlap a requested age range."""
    available = [(20, 29), (30, 39), (40, 49), (50, 59), (60, 69), (70, 79)]
    low = minimum if minimum is not None else 20
    high = maximum if maximum is not None else 79
    return [f"{start}-{end}" for start, end in available if start <= high and end >= low]


class GTExClient:
    """Retrieve public GTEx sample metadata while retaining acquisition context."""

    def __init__(self, client: httpx.Client | None = None) -> None:
        self._client = client or httpx.Client(timeout=60, follow_redirects=True)

    def fetch(
        self, panel: ReferencePanelSpec, items_per_page: int = 1000
    ) -> tuple[list[SampleRecord], ProvenanceRecord]:
        brackets = age_brackets(panel.age.minimum, panel.age.maximum)
        params: dict[str, Any] = {
            "datasetId": GTEX_DATASET_ID,
            "tissueSiteDetailId": panel.tissue,
            "ageBracket": brackets,
            "dataType": "RNASEQ",
            "itemsPerPage": items_per_page,
        }
        rows: list[dict[str, Any]] = []
        page = 0
        while True:
            params["page"] = page
            response = self._client.get(GTEX_SAMPLE_ENDPOINT, params=params)
            response.raise_for_status()
            payload = response.json()
            rows.extend(payload.get("data", []))
            pages = payload.get("paging_info", {}).get("numberOfPages", 1)
            page += 1
            if page >= pages:
                break

        samples = [self._normalize(row, panel) for row in rows]
        provenance = ProvenanceRecord(
            source=DataSource.GTEX,
            endpoint=GTEX_SAMPLE_ENDPOINT,
            retrieved_at=datetime.now(UTC).isoformat(),
            query={key: value for key, value in params.items() if key != "page"},
            source_release=GTEX_DATASET_ID,
            citation_url=GTEX_CITATION_URL,
        )
        return samples, provenance

    @staticmethod
    def _normalize(row: dict[str, Any], panel: ReferencePanelSpec) -> SampleRecord:
        categories = row.get("pathologyNotesCategories") or {}
        present = sorted(name for name, value in categories.items() if value)
        record = SampleRecord(
            source=DataSource.GTEX,
            cohort_role="reference",
            cohort_name=panel.name,
            project_id=row.get("datasetId", GTEX_DATASET_ID),
            case_id=row["subjectId"],
            sample_id=row["sampleId"],
            sample_type="Reference tissue",
            primary_site=row.get("tissueSiteDetail"),
            tissue_or_organ_of_origin=row.get("uberonId"),
            age_bracket=row.get("ageBracket"),
            sex_at_birth=row.get("sex"),
            reference_context=panel.context,
            acquisition_context="postmortem tissue donation",
            assay_type=row.get("dataType"),
            ischemic_time_minutes=row.get("ischemicTime"),
            rin=row.get("rin"),
            hardy_scale=row.get("hardyScale"),
            autolysis_score=row.get("autolysisScore"),
            pathology_notes=row.get("pathologyNotes"),
            pathology_categories_present=present,
            source_payload=row,
        )
        required = ["age_bracket", "sex_at_birth", "primary_site", "assay_type"]
        record.metadata_missing = [name for name in required if getattr(record, name) is None]
        return record
