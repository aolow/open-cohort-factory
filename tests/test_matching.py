from open_cohort_factory.matching import exact_match
from open_cohort_factory.models import DataSource, MatchingSpec, ReferenceContext, SampleRecord


def donor(role: str, donor_id: str, age: str, sex: str) -> SampleRecord:
    source = DataSource.GDC if role == "disease" else DataSource.GTEX
    return SampleRecord(
        source=source,
        cohort_role=role,  # type: ignore[arg-type]
        cohort_name="disease_cohort" if role == "disease" else "GTEx lung",
        project_id="project",
        case_id=donor_id,
        sample_id=f"{donor_id}-sample",
        age_bracket=age,
        sex_at_birth=sex,
    )


def test_exact_match_is_balanced_deterministic_and_auditable() -> None:
    samples = [
        donor("disease", "d1", "60-69", "female"),
        donor("disease", "d2", "60-69", "female"),
        donor("disease", "d3", "80-89", "male"),
        donor("reference", "r1", "60-69", "female"),
    ]
    spec = MatchingSpec(enabled=True, seed=11)
    first = exact_match(samples, spec)[0]
    second = exact_match(samples, spec)[0]
    assert first == second
    assert first["matched_pairs"] == 1
    assert len(first["selected_donors"]) == 2
    assert len(first["excluded_donors"]) == 2
    assert {row["reason"] for row in first["excluded_donors"]} == {
        "no_opposite_donors_in_stratum",
        "surplus_in_stratum",
    }


def test_matching_can_be_disabled() -> None:
    assert exact_match([], MatchingSpec(enabled=False)) == []


def test_adjacent_reference_uses_within_donor_linkage() -> None:
    disease = donor("disease", "p1", "60-69", "female")
    reference = donor("reference", "p1", "60-69", "female").model_copy(
        update={
            "source": DataSource.GDC,
            "cohort_name": "adjacent",
            "reference_context": ReferenceContext.ADJACENT_NON_TUMOR,
        }
    )
    result = exact_match([disease, reference], MatchingSpec(enabled=True))[0]
    assert result["method"] == "within_donor"
    assert result["matched_pairs"] == 1
    assert {row["donor_id"] for row in result["selected_donors"]} == {"p1"}
