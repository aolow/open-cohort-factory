"""Expression summaries that keep unmatched and matched analyses explicit."""

from __future__ import annotations

import math
import random
import statistics
from collections import defaultdict
from typing import Any


def summarize_expression(
    rows: list[dict[str, Any]], matching: list[dict[str, Any]]
) -> dict[str, Any]:
    """Summarize all eligible samples and each donor-matched comparison."""
    disease_rows = [row for row in rows if row["cohort_role"] == "disease"]
    reference_names = sorted(
        {row["cohort_name"] for row in rows if row["cohort_role"] == "reference"}
    )
    analyses = [
        _comparison(
            f"all eligible: {reference_name}",
            disease_rows + [row for row in rows if row["cohort_name"] == reference_name],
        )
        for reference_name in reference_names
    ]
    for reference_name in reference_names:
        reference_rows = [row for row in rows if row["cohort_name"] == reference_name]
        if not any(row.get("reference_context") == "adjacent_non_tumor" for row in reference_rows):
            continue
        disease_donors = {row["donor_id"] for row in disease_rows}
        reference_donors = {row["donor_id"] for row in reference_rows}
        paired_donors = disease_donors & reference_donors
        paired_rows = [
            row for row in disease_rows + reference_rows if row["donor_id"] in paired_donors
        ]
        pairs = [(donor_id, donor_id) for donor_id in sorted(paired_donors)]
        comparison = _comparison(f"within donor: {reference_name}", paired_rows, donor_pairs=pairs)
        comparison["within_donor_pairs"] = len(paired_donors)
        if comparison["genes"]:
            analyses.append(comparison)
    available = {(row["cohort_role"], row["donor_id"]) for row in rows}
    for result in matching:
        reference_rows = [row for row in rows if row["cohort_name"] == result["reference_panel"]]
        if any(row.get("reference_context") == "adjacent_non_tumor" for row in reference_rows):
            continue
        pair_members: defaultdict[str, set[tuple[str, str]]] = defaultdict(set)
        for row in result.get("selected_donors", []):
            pair_members[row["pair_id"]].add((row["cohort_role"], row["donor_id"]))
        complete_pairs = {
            pair_id
            for pair_id, members in pair_members.items()
            if len(members) == 2 and members <= available
        }
        selected = set().union(*(pair_members[pair_id] for pair_id in complete_pairs))
        pairs = []
        for pair_id in sorted(complete_pairs):
            members = pair_members[pair_id]
            disease_id = next(donor for role, donor in members if role == "disease")
            reference_id = next(donor for role, donor in members if role == "reference")
            pairs.append((disease_id, reference_id))
        matched_rows = [
            row
            for row in disease_rows + reference_rows
            if (row["cohort_role"], row["donor_id"]) in selected
        ]
        comparison = _comparison(
            f"matched: {result['reference_panel']}", matched_rows, donor_pairs=pairs
        )
        comparison["complete_expression_pairs"] = len(complete_pairs)
        if comparison["genes"]:
            analyses.append(comparison)
    return {
        "scale": "log2(TPM + 1)",
        "comparisons": analyses,
        "sensitivity_chart": _sensitivity_chart(analyses),
    }


def summarize_pan_tissue(rows: list[dict[str, Any]], tpm_threshold: float) -> dict[str, Any]:
    """Create donor-aware normal-tissue safety summaries for each requested gene."""
    log_threshold = math.log2(tpm_threshold + 1)
    donor_values: defaultdict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in rows:
        donor_values[(row["gene"], row["tissue"], row["donor_id"])].append(row["log2_tpm"])
    grouped: defaultdict[tuple[str, str], list[float]] = defaultdict(list)
    for (gene, tissue, _), values in donor_values.items():
        grouped[(gene, tissue)].append(statistics.median(values))

    gene_summaries = []
    for gene in sorted({gene for gene, _ in grouped}):
        tissues: list[dict[str, Any]] = []
        for row_gene, tissue in sorted(grouped):
            if row_gene != gene:
                continue
            values = grouped[(gene, tissue)]
            tissues.append(
                {
                    "tissue": tissue,
                    "donors": len(values),
                    "median": round(statistics.median(values), 4),
                    "q1": round(_percentile(values, 0.25), 4),
                    "q3": round(_percentile(values, 0.75), 4),
                    "fraction_above_threshold": round(
                        sum(value >= log_threshold for value in values) / len(values), 3
                    ),
                }
            )
        ranked = sorted(tissues, key=lambda item: item["median"], reverse=True)
        maximum = max((item["median"] for item in ranked), default=1) or 1
        top_tissues = [
            {**item, "chart_percent": round(item["median"] / maximum * 100, 2)}
            for item in ranked[:10]
        ]
        gene_summaries.append(
            {
                "gene": gene,
                "tissues": ranked,
                "top_tissues": top_tissues,
                "tissues_with_majority_above_threshold": sum(
                    item["fraction_above_threshold"] >= 0.5 for item in ranked
                ),
            }
        )
    return {
        "scale": "log2(TPM + 1)",
        "tpm_threshold": tpm_threshold,
        "genes": gene_summaries,
    }


def summarize_cell_type_expression(
    rows: list[dict[str, Any]], minimum_donors: int = 3, minimum_cells: int = 50
) -> dict[str, Any]:
    """Rank cell types using donor-level detection and expression summaries."""
    genes = []
    for gene in sorted({str(row["gene"]) for row in rows}):
        gene_rows = [row for row in rows if row["gene"] == gene]
        cell_types = []
        for cell_type in sorted({str(row["cell_type"]) for row in gene_rows}):
            typed = [row for row in gene_rows if row["cell_type"] == cell_type]
            detection = [float(row["fraction_detected"]) for row in typed]
            expression = [float(row["mean_log1p_raw_count"]) for row in typed]
            donor_count = len({row["donor_id"] for row in typed})
            cell_count = sum(int(row["cells"]) for row in typed)
            leave_one_out = [
                statistics.median(
                    float(row["fraction_detected"]) for row in typed if row["donor_id"] != excluded
                )
                for excluded in {row["donor_id"] for row in typed}
                if donor_count > 1
            ]
            cell_types.append(
                {
                    "cell_type": cell_type,
                    "ontology_term_id": typed[0]["cell_type_ontology_term_id"],
                    "donors": donor_count,
                    "cells": cell_count,
                    "median_fraction_detected": round(statistics.median(detection), 3),
                    "median_mean_log1p_raw_count": round(statistics.median(expression), 4),
                    "support": (
                        "supported"
                        if donor_count >= minimum_donors and cell_count >= minimum_cells
                        else "limited"
                    ),
                    "leave_one_donor_out_detection_min": (
                        round(min(leave_one_out), 3) if leave_one_out else None
                    ),
                    "leave_one_donor_out_detection_max": (
                        round(max(leave_one_out), 3) if leave_one_out else None
                    ),
                }
            )
        ranked = sorted(
            cell_types,
            key=lambda item: (
                item["median_fraction_detected"],
                item["median_mean_log1p_raw_count"],
            ),
            reverse=True,
        )
        donors = sorted({row["donor_id"] for row in gene_rows})
        full_top = ranked[0]["cell_type"] if ranked else None
        leave_one_out_tops = [
            _top_cell_type([row for row in gene_rows if row["donor_id"] != donor])
            for donor in donors
            if len(donors) > 1
        ]
        genes.append(
            {
                "gene": gene,
                "cell_types": ranked,
                "top_cell_types": ranked[:8],
                "top_cell_type": full_top,
                "leave_one_donor_out_top_agreement": (
                    round(sum(value == full_top for value in leave_one_out_tops) / len(donors), 3)
                    if leave_one_out_tops
                    else None
                ),
            }
        )
    return {
        "genes": genes,
        "independence_unit": "donor",
        "ranking": "median donor detection fraction, then median donor mean log1p raw count",
        "support_thresholds": {"minimum_donors": minimum_donors, "minimum_cells": minimum_cells},
    }


def _top_cell_type(rows: list[dict[str, Any]]) -> str | None:
    candidates = []
    for cell_type in {str(row["cell_type"]) for row in rows}:
        typed = [row for row in rows if row["cell_type"] == cell_type]
        candidates.append(
            (
                statistics.median(float(row["fraction_detected"]) for row in typed),
                statistics.median(float(row["mean_log1p_raw_count"]) for row in typed),
                cell_type,
            )
        )
    return max(candidates)[2] if candidates else None


def _sensitivity_chart(analyses: list[dict[str, Any]]) -> dict[str, Any]:
    values = [
        max(abs(row["median_difference"]), abs(row["ci_lower"]), abs(row["ci_upper"]))
        for analysis in analyses
        for row in analysis["genes"]
    ]
    limit = max(1, math.ceil(max(values, default=1)))
    genes = sorted({row["gene"] for analysis in analyses for row in analysis["genes"]})
    return {
        "minimum": -limit,
        "maximum": limit,
        "genes": [
            {
                "gene": gene,
                "estimates": [
                    {
                        "comparison": analysis["name"],
                        "value": row["median_difference"],
                        "ci_lower": row["ci_lower"],
                        "ci_upper": row["ci_upper"],
                        "position": round(50 + row["median_difference"] / (2 * limit) * 100, 2),
                        "ci_left": round(50 + row["ci_lower"] / (2 * limit) * 100, 2),
                        "ci_width": round(
                            (row["ci_upper"] - row["ci_lower"]) / (2 * limit) * 100, 2
                        ),
                    }
                    for analysis in analyses
                    for row in analysis["genes"]
                    if row["gene"] == gene
                ],
            }
            for gene in genes
        ],
    }


def _comparison(
    name: str,
    rows: list[dict[str, Any]],
    donor_pairs: list[tuple[str, str]] | None = None,
) -> dict[str, Any]:
    sample_values: defaultdict[tuple[str, str], list[float]] = defaultdict(list)
    donor_values: defaultdict[tuple[str, str, str], list[float]] = defaultdict(list)
    donors: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        key = (row["gene"], row["cohort_role"])
        sample_values[key].append(row["log2_tpm"])
        donor_values[(row["gene"], row["cohort_role"], row["donor_id"])].append(row["log2_tpm"])
        donors[key].add(row["donor_id"])

    genes = sorted({gene for gene, _ in sample_values})
    summaries: list[dict[str, Any]] = []
    for gene in genes:
        tumor_samples = sample_values[(gene, "disease")]
        reference_samples = sample_values[(gene, "reference")]
        tumor_by_donor = {
            donor: statistics.median(values)
            for (row_gene, role, donor), values in donor_values.items()
            if row_gene == gene and role == "disease"
        }
        reference_by_donor = {
            donor: statistics.median(values)
            for (row_gene, role, donor), values in donor_values.items()
            if row_gene == gene and role == "reference"
        }
        tumor = list(tumor_by_donor.values())
        reference = list(reference_by_donor.values())
        if not tumor or not reference:
            continue
        tumor_median = statistics.median(tumor)
        reference_median = statistics.median(reference)
        interval = _bootstrap_interval(
            tumor_by_donor, reference_by_donor, donor_pairs, seed=f"2026|{name}|{gene}"
        )
        summaries.append(
            {
                "gene": gene,
                "disease_samples": len(tumor_samples),
                "reference_samples": len(reference_samples),
                "disease_donors": len(donors[(gene, "disease")]),
                "reference_donors": len(donors[(gene, "reference")]),
                "disease_median": round(tumor_median, 4),
                "reference_median": round(reference_median, 4),
                "disease_q1": round(_percentile(tumor, 0.25), 4),
                "disease_q3": round(_percentile(tumor, 0.75), 4),
                "reference_q1": round(_percentile(reference, 0.25), 4),
                "reference_q3": round(_percentile(reference, 0.75), 4),
                "median_difference": round(tumor_median - reference_median, 4),
                "ci_lower": round(interval[0], 4),
                "ci_upper": round(interval[1], 4),
            }
        )
    return {"name": name, "genes": summaries}


def _percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def _bootstrap_interval(
    disease: dict[str, float],
    reference: dict[str, float],
    donor_pairs: list[tuple[str, str]] | None,
    seed: str,
    iterations: int = 1000,
) -> tuple[float, float]:
    rng = random.Random(seed)
    estimates: list[float] = []
    complete = [
        (disease[left], reference[right])
        for left, right in (donor_pairs or [])
        if left in disease and right in reference
    ]
    if complete:
        for _ in range(iterations):
            resampled = rng.choices(complete, k=len(complete))
            estimates.append(
                statistics.median(left for left, _ in resampled)
                - statistics.median(right for _, right in resampled)
            )
    else:
        disease_values = list(disease.values())
        reference_values = list(reference.values())
        for _ in range(iterations):
            estimates.append(
                statistics.median(rng.choices(disease_values, k=len(disease_values)))
                - statistics.median(rng.choices(reference_values, k=len(reference_values)))
            )
    return _percentile(estimates, 0.025), _percentile(estimates, 0.975)
