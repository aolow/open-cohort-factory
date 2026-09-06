import httpx

from open_cohort_factory.models import (
    AgeRange,
    DataSource,
    DiseaseCohortSpec,
    ReferenceContext,
    ReferencePanelSpec,
)
from open_cohort_factory.sources.gdc import GDCClient, build_filters


def test_build_filters_includes_age_and_sample_type() -> None:
    spec = DiseaseCohortSpec(
        projects=["TCGA-LUAD"], age=AgeRange(minimum=60), sample_types=["Primary Tumor"]
    )
    filters = build_filters(spec)
    fields = [clause["content"]["field"] for clause in filters["content"]]
    assert "demographic.age_at_index" in fields
    assert "samples.sample_type" in fields


def test_normalizes_samples_without_calling_adjacent_tissue_healthy() -> None:
    payload = {
        "data": {
            "hits": [
                {
                    "case_id": "case-1",
                    "primary_site": "Lung",
                    "project": {"project_id": "TCGA-LUAD"},
                    "demographic": {"age_at_index": 67, "sex_at_birth": "female"},
                    "diagnoses": [{"tissue_or_organ_of_origin": "Upper lobe, lung"}],
                    "samples": [{"sample_id": "sample-1", "sample_type": "Solid Tissue Normal"}],
                }
            ]
        }
    }

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=payload, request=request)

    client = GDCClient(httpx.Client(transport=httpx.MockTransport(handler)))
    spec = DiseaseCohortSpec(projects=["TCGA-LUAD"], sample_types=["Solid Tissue Normal"])
    samples, _ = client.fetch(spec)
    assert samples[0].reference_context == ReferenceContext.ADJACENT_NON_TUMOR
    assert samples[0].acquisition_context == "cancer study participant"


def test_gdc_reference_is_not_normalized_as_disease() -> None:
    payload = {
        "data": {
            "hits": [{
                "case_id": "case-1",
                "primary_site": "Bronchus and lung",
                "project": {"project_id": "TCGA-LUAD"},
                "demographic": {"age_at_index": 67, "sex_at_birth": "female"},
                "samples": [{"sample_id": "normal-1", "sample_type": "Solid Tissue Normal"}],
            }]
        }
    }
    client = GDCClient(httpx.Client(transport=httpx.MockTransport(
        lambda request: httpx.Response(200, json=payload, request=request)
    )))
    panel = ReferencePanelSpec(
        name="adjacent",
        source=DataSource.GDC,
        tissue="Bronchus and lung",
        context=ReferenceContext.ADJACENT_NON_TUMOR,
        sample_types=["Solid Tissue Normal"],
    )
    samples, provenance = client.fetch_reference(
        panel, DiseaseCohortSpec(projects=["TCGA-LUAD"])
    )
    assert samples[0].cohort_role == "reference"
    assert samples[0].cohort_name == "adjacent"
    assert provenance.query["reference_panel"] == "adjacent"
