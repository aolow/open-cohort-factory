from open_cohort_factory.audit import methodological_warnings, summarize
from open_cohort_factory.models import (
    DataSource,
    DiseaseCohortSpec,
    ProjectSpec,
    SampleRecord,
)


def test_summary_counts_donors_not_just_samples() -> None:
    spec = ProjectSpec(
        name="test", description="test", disease_cohort=DiseaseCohortSpec(projects=["P"])
    )
    samples = [
        SampleRecord(
            source=DataSource.GDC,
            cohort_role="disease",
            cohort_name="disease_cohort",
            project_id="P",
            case_id="C",
            sample_id="S1",
        ),
        SampleRecord(
            source=DataSource.GDC,
            cohort_role="disease",
            cohort_name="disease_cohort",
            project_id="P",
            case_id="C",
            sample_id="S2",
        ),
    ]
    result = summarize(samples, spec)
    assert result["sample_count"] == 2
    assert result["donor_count"] == 1
    assert "universal normal baseline" in methodological_warnings(spec, samples)[0]


def test_covariate_landscape_flags_source_confounded_availability() -> None:
    spec = ProjectSpec(
        name="test", description="test", disease_cohort=DiseaseCohortSpec(projects=["P"])
    )
    samples = [
        SampleRecord(
            source=DataSource.GDC,
            cohort_role="disease",
            cohort_name="disease",
            project_id="P",
            case_id="d1",
            sample_id="d1-s",
            race="white",
        ),
        SampleRecord(
            source=DataSource.GTEX,
            cohort_role="reference",
            cohort_name="GTEx",
            project_id="G",
            case_id="r1",
            sample_id="r1-s",
        ),
    ]
    fields = summarize(samples, spec)["covariate_landscape"]["fields"]
    race = next(field for field in fields if field["field"] == "race")
    assert race["availability_status"] == "source_confounded"
