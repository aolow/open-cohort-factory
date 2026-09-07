from open_cohort_factory.models import DataSource, ReferenceContext, SampleRecord
from open_cohort_factory.reference_fitness import assess_reference_fitness


def test_reference_fitness_preserves_dimensions_instead_of_composite_score() -> None:
    sample = SampleRecord(
        source=DataSource.GTEX,
        cohort_role="reference",
        cohort_name="GTEx lung",
        project_id="gtex",
        case_id="r1",
        sample_id="s1",
        reference_context=ReferenceContext.POSTMORTEM_REFERENCE,
    )
    coverage = [
        {"field": field, "coverage": [{"population": "GTEx lung", "percent": percent}]}
        for field, percent in [
            ("age_bracket", 100),
            ("sex_at_birth", 100),
            ("race", 0),
            ("ethnicity", 0),
        ]
    ]
    summary = {
        "populations": [{"role": "reference", "name": "GTEx lung", "donors": 1}],
        "comparability": [
            {
                "reference_panel": "GTEx lung",
                "uncovered_disease_age_brackets": ["80-89"],
                "female_proportion": {"absolute_difference": 0.2},
                "decision": "insufficient_reference_donors",
            }
        ],
        "covariate_landscape": {"fields": coverage},
    }
    result = assess_reference_fitness([sample], summary, [], None)[0]
    assert "score" not in result
    statuses = {row["dimension"]: row["status"] for row in result["dimensions"]}
    assert statuses["Population support"] == "limited"
    assert statuses["Baseline metadata"] == "limited"
    assert "unqualified healthy counterfactual" in result["not_supported"]
