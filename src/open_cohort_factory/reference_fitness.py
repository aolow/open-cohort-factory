"""Multi-dimensional reference-panel fitness assessments without composite scores."""

from __future__ import annotations

from typing import Any

from .models import SampleRecord


def assess_reference_fitness(
    samples: list[SampleRecord],
    summary: dict[str, Any],
    matching: list[dict[str, Any]],
    expression: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Describe panel suitability by dimension and intended analytical use."""
    comparison_by_name = {
        row["reference_panel"]: row for row in summary.get("comparability", [])
    }
    matching_by_name = {row["reference_panel"]: row for row in matching}
    population_by_name = {
        row["name"]: row for row in summary.get("populations", []) if row["role"] == "reference"
    }
    expression_by_name = _expression_coverage(expression)
    results: list[dict[str, Any]] = []

    reference_names = sorted(
        {sample.cohort_name for sample in samples if sample.cohort_role == "reference"}
    )
    for name in reference_names:
        members = [sample for sample in samples if sample.cohort_name == name]
        context = members[0].reference_context.value if members[0].reference_context else "unknown"
        population = population_by_name[name]
        comparison = comparison_by_name[name]
        match = matching_by_name.get(name)
        expression_coverage = expression_by_name.get(name)
        dimensions = [
            _population_dimension(population["donors"]),
            _demographic_dimension(comparison),
            _acquisition_dimension(context),
            _matching_dimension(match, comparison),
            _measurement_dimension(expression_coverage, population["donors"], members),
            _metadata_dimension(name, summary["covariate_landscape"]),
        ]
        results.append(
            {
                "reference_panel": name,
                "context": context,
                "dimensions": dimensions,
                "recommended_uses": _recommended_uses(context, dimensions, members),
                "not_supported": _unsupported_uses(context, dimensions),
            }
        )
    return results


def _population_dimension(donors: int) -> dict[str, Any]:
    if donors >= 50:
        status = "supported"
    elif donors >= 10:
        status = "caution"
    else:
        status = "limited"
    return {
        "dimension": "Population support",
        "status": status,
        "evidence": f"{donors} eligible donors",
    }


def _demographic_dimension(comparison: dict[str, Any]) -> dict[str, Any]:
    uncovered = comparison["uncovered_disease_age_brackets"]
    sex_difference = comparison["female_proportion"]["absolute_difference"]
    status = "supported"
    reasons = []
    if uncovered:
        status = "limited"
        reasons.append(f"uncovered disease ages: {', '.join(uncovered)}")
    if sex_difference is None:
        status = "not_evaluable"
        reasons.append("sex distribution unavailable")
    elif sex_difference >= 0.1:
        status = "caution" if status == "supported" else status
        reasons.append(f"female-proportion difference {sex_difference}")
    return {
        "dimension": "Demographic alignment",
        "status": status,
        "evidence": "; ".join(reasons) or "no large observed age/sex mismatch",
    }


def _acquisition_dimension(context: str) -> dict[str, Any]:
    mapping = {
        "adjacent_non_tumor": (
            "caution",
            "same disease setting, but tissue may carry field or systemic cancer effects",
        ),
        "postmortem_reference": (
            "caution",
            "postmortem procurement differs from surgical tumor collection",
        ),
        "organ_donor": (
            "caution",
            "organ-donor procurement and eligibility differ from the disease cohort",
        ),
    }
    status, evidence = mapping.get(context, ("not_evaluable", "acquisition alignment unknown"))
    return {"dimension": "Acquisition alignment", "status": status, "evidence": evidence}


def _matching_dimension(
    match: dict[str, Any] | None, comparison: dict[str, Any]
) -> dict[str, Any]:
    if not match:
        return {
            "dimension": "Design linkage",
            "status": "not_evaluable",
            "evidence": "no matching analysis",
        }
    pairs = match["matched_pairs"]
    if match["method"] == "within_donor":
        evidence = f"{pairs} within-participant tumor/reference pairs"
        status = "supported" if pairs >= 10 else "limited"
    else:
        evidence = f"{pairs} exact demographic pairs; residual covariates remain"
        status = "supported" if pairs >= 50 else "caution"
    if comparison["decision"] == "insufficient_reference_donors":
        status = "limited"
    return {"dimension": "Design linkage", "status": status, "evidence": evidence}


def _measurement_dimension(
    coverage: int | None, donors: int, members: list[SampleRecord]
) -> dict[str, Any]:
    if any(sample.cell_type is not None for sample in members):
        return {
            "dimension": "Measurement overlap",
            "status": "not_evaluable",
            "evidence": "single-cell context is not on the bulk Xena measurement scale",
        }
    if coverage is None:
        return {
            "dimension": "Measurement overlap",
            "status": "not_evaluable",
            "evidence": "no expression measurements returned",
        }
    fraction = coverage / donors if donors else 0
    status = "supported" if fraction >= 0.8 else "caution" if fraction >= 0.5 else "limited"
    return {
        "dimension": "Measurement overlap",
        "status": status,
        "evidence": (
            f"{coverage}/{donors} donors represented in the expression slice ({fraction:.0%})"
        ),
    }


def _metadata_dimension(name: str, landscape: dict[str, Any]) -> dict[str, Any]:
    baseline = {"age_bracket", "sex_at_birth", "race", "ethnicity"}
    rows = [row for row in landscape["fields"] if row["field"] in baseline]
    coverage = {
        row["field"]: next(
            item["percent"] for item in row["coverage"] if item["population"] == name
        )
        for row in rows
    }
    low = [field for field, percent in coverage.items() if percent < 50]
    status = "supported" if not low else "limited"
    return {
        "dimension": "Baseline metadata",
        "status": status,
        "evidence": (
            "age, sex, race, and ethnicity mostly observed"
            if not low
            else "<50% observed: " + ", ".join(low)
        ),
    }


def _expression_coverage(expression: dict[str, Any] | None) -> dict[str, int]:
    if not expression:
        return {}
    coverage: dict[str, int] = {}
    for analysis in expression["comparisons"]:
        prefix = "all eligible: "
        if analysis["name"].startswith(prefix) and analysis["genes"]:
            coverage[analysis["name"][len(prefix) :]] = analysis["genes"][0]["reference_donors"]
    return coverage


def _recommended_uses(
    context: str, dimensions: list[dict[str, Any]], members: list[SampleRecord]
) -> list[str]:
    uses = ["reference-context sensitivity analysis"]
    if context == "adjacent_non_tumor":
        uses.append("within-participant tumor-versus-adjacent description")
    if context == "postmortem_reference":
        uses.append("population-level adult tissue expression context")
    if any(sample.cell_type is not None for sample in members):
        uses.append("cell-type attribution and hypothesis generation")
    if all(row["status"] != "limited" for row in dimensions):
        uses.append("primary descriptive comparator with stated caveats")
    return uses


def _unsupported_uses(context: str, dimensions: list[dict[str, Any]]) -> list[str]:
    uses = ["unqualified healthy counterfactual", "causal disease-effect estimate"]
    if any(row["status"] == "limited" for row in dimensions):
        uses.append("sole primary reference without sensitivity analysis")
    if context == "organ_donor":
        uses.append("bulk-expression comparator")
    return uses
