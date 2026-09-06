from open_cohort_factory.comparability import compare_populations
from open_cohort_factory.models import DataSource, ReferenceContext, SampleRecord


def record(
    source: DataSource,
    role: str,
    name: str,
    donor: str,
    age_bracket: str,
    sex: str,
) -> SampleRecord:
    return SampleRecord(
        source=source,
        cohort_role=role,  # type: ignore[arg-type]
        cohort_name=name,
        project_id="project",
        case_id=donor,
        sample_id=f"{donor}-sample",
        age_bracket=age_bracket,
        sex_at_birth=sex,
        reference_context=(ReferenceContext.POSTMORTEM_REFERENCE if role == "reference" else None),
    )


def test_flags_age_coverage_and_sex_imbalance_at_donor_level() -> None:
    samples = [
        record(DataSource.GDC, "disease", "disease_cohort", "d1", "60-69", "female"),
        record(DataSource.GDC, "disease", "disease_cohort", "d2", "80-89", "female"),
        record(DataSource.GTEX, "reference", "GTEx lung", "r1", "60-69", "male"),
        record(DataSource.GTEX, "reference", "GTEx lung", "r2", "60-69", "male"),
    ]
    comparison = compare_populations(samples)[0]
    assert comparison["uncovered_disease_age_brackets"] == ["80-89"]
    assert comparison["female_proportion"]["absolute_difference"] == 1.0
    assert comparison["decision"] == "context_only_for_full_cohort"
