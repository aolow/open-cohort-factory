"""Canonical, source-independent cohort models."""

from __future__ import annotations

from enum import StrEnum
from typing import Any, Literal, TypeAlias

from pydantic import BaseModel, Field, model_validator


class ReferenceContext(StrEnum):
    """How tissue entered a reference panel.

    The values intentionally avoid using ``normal`` as an unqualified biological state.
    """

    ADJACENT_NON_TUMOR = "adjacent_non_tumor"
    POSTMORTEM_REFERENCE = "postmortem_reference"
    ORGAN_DONOR = "organ_donor"
    SURGICAL_NON_DISEASED = "surgical_non_diseased"
    HEALTHY_VOLUNTEER = "healthy_volunteer"
    DISEASE_CONTROL = "disease_control"
    UNKNOWN = "unknown_or_incompletely_characterized"


class DataSource(StrEnum):
    GDC = "gdc"
    GTEX = "gtex"
    TABULA_SAPIENS = "tabula_sapiens"
    CELLXGENE = "cellxgene"


class AgeRange(BaseModel):
    minimum: int | None = Field(default=None, ge=0, le=120)
    maximum: int | None = Field(default=None, ge=0, le=120)

    @model_validator(mode="after")
    def ordered(self) -> AgeRange:
        if self.minimum is not None and self.maximum is not None and self.minimum > self.maximum:
            raise ValueError("minimum age cannot exceed maximum age")
        return self


class DiseaseCohortSpec(BaseModel):
    source: Literal[DataSource.GDC] = DataSource.GDC
    projects: list[str] = Field(min_length=1)
    primary_sites: list[str] = Field(default_factory=list)
    age: AgeRange = Field(default_factory=AgeRange)
    sex_at_birth: list[str] = Field(default_factory=list)
    sample_types: list[str] = Field(default_factory=lambda: ["Primary Tumor"])


class ReferencePanelSpec(BaseModel):
    name: str
    source: DataSource
    tissue: str
    context: ReferenceContext
    age: AgeRange = Field(default_factory=AgeRange)
    resolution: Literal["bulk", "cell_type", "single_cell"] = "bulk"
    notes: str | None = None


StratificationField: TypeAlias = Literal["age_group", "sex_at_birth", "source"]


def _default_stratification() -> list[StratificationField]:
    return ["age_group", "sex_at_birth", "source"]


class ComparabilitySpec(BaseModel):
    stratify_by: list[StratificationField] = Field(default_factory=_default_stratification)
    minimum_donors_per_group: int = Field(default=5, ge=1)
    preserve_source_effects: bool = True
    missing_metadata: Literal["report", "exclude"] = "report"


class ProjectSpec(BaseModel):
    name: str
    description: str
    disease_cohort: DiseaseCohortSpec
    reference_panels: list[ReferencePanelSpec] = Field(default_factory=list)
    comparability: ComparabilitySpec = Field(default_factory=ComparabilitySpec)


class SampleRecord(BaseModel):
    source: DataSource
    cohort_role: Literal["disease", "reference"]
    cohort_name: str
    project_id: str
    case_id: str
    sample_id: str
    sample_type: str | None = None
    primary_site: str | None = None
    tissue_or_organ_of_origin: str | None = None
    age_at_index: int | None = None
    age_bracket: str | None = None
    age_at_diagnosis_days: int | None = None
    sex_at_birth: str | None = None
    race: str | None = None
    ethnicity: str | None = None
    vital_status: str | None = None
    reference_context: ReferenceContext | None = None
    acquisition_context: str | None = None
    assay_type: str | None = None
    ischemic_time_minutes: int | None = None
    rin: float | None = None
    hardy_scale: str | None = None
    autolysis_score: str | None = None
    pathology_notes: str | None = None
    pathology_categories_present: list[str] = Field(default_factory=list)
    metadata_missing: list[str] = Field(default_factory=list)
    source_payload: dict[str, Any] = Field(default_factory=dict, exclude=True)


class ProvenanceRecord(BaseModel):
    source: DataSource
    endpoint: str
    retrieved_at: str
    query: dict[str, Any]
    source_release: str | None = None
    citation_url: str
    license_url: str | None = None


class BuildResult(BaseModel):
    specification: ProjectSpec
    samples: list[SampleRecord]
    provenance: list[ProvenanceRecord]
    warnings: list[str] = Field(default_factory=list)
