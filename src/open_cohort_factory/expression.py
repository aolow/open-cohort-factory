"""Expression summaries that keep unmatched and matched analyses explicit."""

from __future__ import annotations

import statistics
from collections import defaultdict
from typing import Any


def summarize_expression(
    rows: list[dict[str, Any]], matching: list[dict[str, Any]]
) -> dict[str, Any]:
    """Summarize all eligible samples and each donor-matched comparison."""
    analyses = [_comparison("all_eligible", rows)]
    available = {(row["cohort_role"], row["donor_id"]) for row in rows}
    for result in matching:
        pair_members: defaultdict[str, set[tuple[str, str]]] = defaultdict(set)
        for row in result.get("selected_donors", []):
            pair_members[row["pair_id"]].add((row["cohort_role"], row["donor_id"]))
        complete_pairs = {
            pair_id
            for pair_id, members in pair_members.items()
            if len(members) == 2 and members <= available
        }
        selected = set().union(*(pair_members[pair_id] for pair_id in complete_pairs))
        matched_rows = [
            row for row in rows if (row["cohort_role"], row["donor_id"]) in selected
        ]
        comparison = _comparison(f"matched: {result['reference_panel']}", matched_rows)
        comparison["complete_expression_pairs"] = len(complete_pairs)
        if comparison["genes"]:
            analyses.append(comparison)
    return {"scale": "log2(TPM + 1)", "comparisons": analyses}


def _comparison(name: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    sample_values: defaultdict[tuple[str, str], list[float]] = defaultdict(list)
    donor_values: defaultdict[tuple[str, str, str], list[float]] = defaultdict(list)
    donors: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        key = (row["gene"], row["cohort_role"])
        sample_values[key].append(row["log2_tpm"])
        donor_values[(row["gene"], row["cohort_role"], row["donor_id"])].append(
            row["log2_tpm"]
        )
        donors[key].add(row["donor_id"])

    genes = sorted({gene for gene, _ in sample_values})
    summaries: list[dict[str, Any]] = []
    for gene in genes:
        tumor_samples = sample_values[(gene, "disease")]
        reference_samples = sample_values[(gene, "reference")]
        tumor = [
            statistics.median(values)
            for (row_gene, role, _), values in donor_values.items()
            if row_gene == gene and role == "disease"
        ]
        reference = [
            statistics.median(values)
            for (row_gene, role, _), values in donor_values.items()
            if row_gene == gene and role == "reference"
        ]
        if not tumor or not reference:
            continue
        tumor_median = statistics.median(tumor)
        reference_median = statistics.median(reference)
        summaries.append(
            {
                "gene": gene,
                "disease_samples": len(tumor_samples),
                "reference_samples": len(reference_samples),
                "disease_donors": len(donors[(gene, "disease")]),
                "reference_donors": len(donors[(gene, "reference")]),
                "disease_median": round(tumor_median, 4),
                "reference_median": round(reference_median, 4),
                "median_difference": round(tumor_median - reference_median, 4),
            }
        )
    return {"name": name, "genes": summaries}
