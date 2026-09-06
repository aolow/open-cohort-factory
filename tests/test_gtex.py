import httpx

from open_cohort_factory.models import DataSource, ReferenceContext, ReferencePanelSpec
from open_cohort_factory.sources.gtex import GTExClient, age_brackets


def test_age_brackets_overlap_requested_range() -> None:
    assert age_brackets(60, 70) == ["60-69", "70-79"]
    assert age_brackets(45, 55) == ["40-49", "50-59"]


def test_gtex_normalization_preserves_pathology_and_acquisition_context() -> None:
    payload = {
        "data": [
            {
                "sampleId": "GTEX-AAAA-0001-SM-ABCDE",
                "subjectId": "GTEX-AAAA",
                "datasetId": "gtex_v10",
                "tissueSiteDetail": "Lung",
                "uberonId": "UBERON:0008952",
                "ageBracket": "60-69",
                "sex": "female",
                "dataType": "RNASEQ",
                "ischemicTime": 480,
                "rin": 7.1,
                "hardyScale": "Fast death - violent",
                "autolysisScore": "Mild",
                "pathologyNotes": "congestion",
                "pathologyNotesCategories": {"congestion": True, "fibrosis": False},
            }
        ],
        "paging_info": {"numberOfPages": 1, "page": 0},
    }

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=payload, request=request)

    panel = ReferencePanelSpec(
        name="GTEx lung",
        source=DataSource.GTEX,
        tissue="Lung",
        context=ReferenceContext.POSTMORTEM_REFERENCE,
        age={"minimum": 60, "maximum": 70},
    )
    client = GTExClient(httpx.Client(transport=httpx.MockTransport(handler)))
    samples, provenance = client.fetch(panel)
    assert samples[0].cohort_role == "reference"
    assert samples[0].pathology_categories_present == ["congestion"]
    assert samples[0].reference_context == ReferenceContext.POSTMORTEM_REFERENCE
    assert provenance.source_release == "gtex_v10"
