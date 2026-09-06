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
