import json

import httpx
import pytest

from open_cohort_factory.models import (
    DataSource,
    ExpressionSpec,
    ReferenceContext,
    SampleRecord,
)
from open_cohort_factory.sources.xena import XenaClient, xena_sample_id


def _sample(source: DataSource, role: str, sample_id: str, submitter_id: str | None = None):
    return SampleRecord(
        source=source,
        cohort_role=role,
        cohort_name="disease" if role == "disease" else "GTEx lung",
        project_id="TCGA-LUAD" if role == "disease" else "gtex_v10",
        case_id=f"donor-{sample_id}",
        sample_id=sample_id,
        sample_submitter_id=submitter_id,
        reference_context=(
            None if role == "disease" else ReferenceContext.POSTMORTEM_REFERENCE
        ),
    )


def test_xena_ids_are_explicit_and_uuid_is_not_guessed():
    gdc = _sample(DataSource.GDC, "disease", "uuid", "TCGA-AA-0001-01A")
    assert xena_sample_id(gdc) == "TCGA-AA-0001-01"
    assert xena_sample_id(_sample(DataSource.GDC, "disease", "uuid")) is None
    assert xena_sample_id(_sample(DataSource.GTEX, "reference", "GTEX-ABC-SM-1")) == (
        "GTEX-ABC-SM-1"
    )


def test_fetch_separates_source_roles_and_records_provenance():
    samples = [
        _sample(DataSource.GDC, "disease", "uuid", "TCGA-AA-0001-01A"),
        _sample(DataSource.GTEX, "reference", "GTEX-ABC-SM-1"),
    ]

    def handler(request: httpx.Request) -> httpx.Response:
        body = request.content.decode()
        assert "TcgaTargetGtex_rsem_gene_tpm" in body
        assert "TCGA-AA-0001-01" in body
        return httpx.Response(
            200,
            json=[{"gene": "EPCAM", "position": [], "scores": [[10.0, 1.0]]}],
        )

    client = httpx.Client(transport=httpx.MockTransport(handler))
    rows, provenance = XenaClient(client=client).fetch(samples, ExpressionSpec(genes=["epcam"]))
    assert [row["source"] for row in rows] == ["gdc", "gtex"]
    assert rows[0]["log2_tpm"] == pytest.approx(10.0014, abs=0.001)
    assert provenance.source == DataSource.XENA_TOIL
    assert provenance.query["input_scale"] == "log2(TPM + 0.001)"
    json.dumps(provenance.model_dump(mode="json"))


def test_fetch_fails_when_only_non_bulk_reference_samples_exist():
    tabula = _sample(DataSource.TABULA_SAPIENS, "reference", "cell-record")
    with pytest.raises(ValueError, match="public Xena identifiers"):
        XenaClient(client=httpx.Client(transport=httpx.MockTransport(lambda request: None))).fetch(
            [tabula], ExpressionSpec(genes=["EPCAM"])
        )
