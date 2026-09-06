from open_cohort_factory.models import DataSource, ReferenceContext, ReferencePanelSpec
from open_cohort_factory.sources.tabula_sapiens import (
    TabulaSapiensClient,
    age_bracket,
    parse_age,
)


def test_age_parsing_and_bracketing() -> None:
    assert parse_age("61-year-old stage") == 61
    assert parse_age("adult stage") is None
    assert age_bracket(61) == "60-69"


def test_normalizes_donor_cell_type_aggregate() -> None:
    panel = ReferencePanelSpec(
        name="Tabula Sapiens lung",
        source=DataSource.TABULA_SAPIENS,
        tissue="lung",
        context=ReferenceContext.ORGAN_DONOR,
        resolution="cell_type",
    )
    row = {
        "dataset_id": "dataset-1",
        "donor_id": "TSP2",
        "sex": "female",
        "age": 61,
        "development_stage": "61-year-old stage",
        "cell_type": "macrophage",
        "cell_type_ontology_term_id": "CL:0000235",
        "tissue": "lung",
        "tissue_ontology_term_id": "UBERON:0002048",
        "disease": "normal",
        "assay": "10x 3' v3",
        "cell_count": 428,
    }
    record = TabulaSapiensClient._normalize(row, panel)
    assert record.case_id == "TSP2"
    assert record.cell_count == 428
    assert record.cell_type == "macrophage"
    assert record.reference_context == ReferenceContext.ORGAN_DONOR
