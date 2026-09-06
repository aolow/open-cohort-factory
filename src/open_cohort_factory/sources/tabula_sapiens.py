"""Donor-aware Tabula Sapiens metadata through the CELLxGENE Census."""

from __future__ import annotations

import importlib
import re
from datetime import UTC, datetime
from typing import Any

from ..models import DataSource, ProvenanceRecord, ReferencePanelSpec, SampleRecord

CENSUS_VERSION = "2025-11-08"
CENSUS_DOCS_URL = "https://chanzuckerberg.github.io/cellxgene-census/"
TABULA_COLLECTION = "Tabula Sapiens"


def parse_age(stage: str) -> int | None:
    match = re.search(r"(\d+)-year-old", stage)
    return int(match.group(1)) if match else None


def age_bracket(age: int) -> str:
    lower = age // 10 * 10
    return f"{lower}-{lower + 9}"


def _matches_age(age: int | None, panel: ReferencePanelSpec) -> bool:
    if age is None:
        return False
    if panel.age.minimum is not None and age < panel.age.minimum:
        return False
    return panel.age.maximum is None or age <= panel.age.maximum


class TabulaSapiensClient:
    """Query cell metadata and emit one record per donor and cell type."""

    def __init__(self, census_version: str = CENSUS_VERSION) -> None:
        self.census_version = census_version

    def fetch(self, panel: ReferencePanelSpec) -> tuple[list[SampleRecord], ProvenanceRecord]:
        try:
            census_api = importlib.import_module("cellxgene_census")
        except ImportError as error:
            raise RuntimeError(
                "Tabula Sapiens requires: uv sync --extra dev --extra single-cell"
            ) from error

        with census_api.open_soma(census_version=self.census_version) as census:
            datasets = (
                census["census_info"]["datasets"]
                .read(
                    column_names=[
                        "dataset_id",
                        "collection_name",
                        "dataset_title",
                        "dataset_total_cell_count",
                    ]
                )
                .concat()
                .to_pandas()
            )
            dataset = self._select_dataset(datasets, panel.tissue)
            dataset_id = str(dataset["dataset_id"])
            columns = [
                "dataset_id",
                "donor_id",
                "sex",
                "development_stage",
                "cell_type",
                "cell_type_ontology_term_id",
                "tissue",
                "tissue_ontology_term_id",
                "disease",
                "assay",
            ]
            cells = census_api.get_obs(
                census,
                "homo_sapiens",
                value_filter=f"dataset_id == '{dataset_id}'",
                column_names=columns,
            )

        cells = cells.copy()
        cells["age"] = cells["development_stage"].map(parse_age)
        cells = cells[cells["age"].map(lambda value: _matches_age(value, panel))]
        group_columns = [
            "dataset_id",
            "donor_id",
            "sex",
            "age",
            "development_stage",
            "cell_type",
            "cell_type_ontology_term_id",
            "tissue",
            "tissue_ontology_term_id",
            "disease",
            "assay",
        ]
        aggregates = (
            cells.groupby(group_columns, observed=True, dropna=False)
            .size()
            .reset_index(name="cell_count")
        )
        samples = [self._normalize(row, panel) for row in aggregates.to_dict(orient="records")]
        provenance = ProvenanceRecord(
            source=DataSource.TABULA_SAPIENS,
            endpoint="CELLxGENE Census soma://census_data/homo_sapiens/obs",
            retrieved_at=datetime.now(UTC).isoformat(),
            query={
                "census_version": self.census_version,
                "collection_name": TABULA_COLLECTION,
                "dataset_id": dataset_id,
                "tissue": panel.tissue,
                "minimum_age": panel.age.minimum,
                "maximum_age": panel.age.maximum,
                "aggregation": "donor_x_cell_type",
            },
            source_release=self.census_version,
            citation_url=CENSUS_DOCS_URL,
        )
        return samples, provenance

    @staticmethod
    def _select_dataset(datasets: Any, tissue: str) -> Any:
        expected = f"tabula sapiens - {tissue.replace('_', ' ')}".lower()
        matches = datasets[
            (datasets["collection_name"].str.lower() == TABULA_COLLECTION.lower())
            & (datasets["dataset_title"].str.replace("_", " ").str.lower() == expected)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one tissue-specific Tabula Sapiens dataset for {tissue!r}; "
                f"found {len(matches)}"
            )
        return matches.iloc[0]

    @staticmethod
    def _normalize(row: dict[str, Any], panel: ReferencePanelSpec) -> SampleRecord:
        age = int(row["age"])
        donor_id = str(row["donor_id"])
        cell_type_id = str(row["cell_type_ontology_term_id"])
        return SampleRecord(
            source=DataSource.TABULA_SAPIENS,
            cohort_role="reference",
            cohort_name=panel.name,
            project_id=str(row["dataset_id"]),
            case_id=donor_id,
            sample_id=f"{row['dataset_id']}:{donor_id}:{cell_type_id}",
            sample_type="Donor cell-type aggregate",
            primary_site=str(row["tissue"]),
            tissue_or_organ_of_origin=str(row["tissue_ontology_term_id"]),
            age_at_index=age,
            age_bracket=age_bracket(age),
            sex_at_birth=str(row["sex"]),
            reference_context=panel.context,
            acquisition_context="organ donor tissue",
            assay_type=str(row["assay"]),
            cell_type=str(row["cell_type"]),
            cell_type_ontology_term_id=cell_type_id,
            cell_count=int(row["cell_count"]),
            development_stage=str(row["development_stage"]),
            disease_label=str(row["disease"]),
        )
