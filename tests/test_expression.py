from open_cohort_factory.expression import (
    summarize_cell_type_expression,
    summarize_expression,
    summarize_pan_tissue,
)


def test_expression_reports_unmatched_and_donor_matched_results():
    rows = [
        dict(
            gene="EPCAM", cohort_role="disease", cohort_name="disease", donor_id="d1", log2_tpm=8.0
        ),
        dict(
            gene="EPCAM", cohort_role="disease", cohort_name="disease", donor_id="d2", log2_tpm=4.0
        ),
        dict(
            gene="EPCAM",
            cohort_role="reference",
            cohort_name="GTEx lung",
            donor_id="r1",
            log2_tpm=2.0,
        ),
        dict(
            gene="EPCAM",
            cohort_role="reference",
            cohort_name="GTEx lung",
            donor_id="r2",
            log2_tpm=0.0,
        ),
    ]
    matching = [
        {
            "reference_panel": "GTEx lung",
            "selected_donors": [
                {"pair_id": "pair-1", "cohort_role": "disease", "donor_id": "d1"},
                {"pair_id": "pair-1", "cohort_role": "reference", "donor_id": "r1"},
            ],
        }
    ]
    result = summarize_expression(rows, matching)
    assert result["comparisons"][0]["genes"][0]["median_difference"] == 5.0
    assert result["comparisons"][1]["genes"][0]["median_difference"] == 6.0
    assert result["comparisons"][1]["complete_expression_pairs"] == 1
    assert result["comparisons"][0]["genes"][0]["ci_lower"] <= 5.0
    assert result["comparisons"][0]["genes"][0]["ci_upper"] >= 5.0
    positions = result["sensitivity_chart"]["genes"][0]["estimates"]
    assert all(0 <= estimate["position"] <= 100 for estimate in positions)


def test_adjacent_reference_adds_within_donor_comparison():
    rows = [
        dict(gene="G", cohort_role="disease", cohort_name="disease", donor_id="p1", log2_tpm=5.0),
        dict(
            gene="G",
            cohort_role="reference",
            cohort_name="adjacent",
            donor_id="p1",
            log2_tpm=2.0,
            reference_context="adjacent_non_tumor",
        ),
        dict(
            gene="G",
            cohort_role="reference",
            cohort_name="adjacent",
            donor_id="p2",
            log2_tpm=1.0,
            reference_context="adjacent_non_tumor",
        ),
    ]
    result = summarize_expression(rows, [])
    within = result["comparisons"][1]
    assert within["within_donor_pairs"] == 1
    assert within["genes"][0]["median_difference"] == 3.0


def test_pan_tissue_summary_is_donor_aware_and_ranked():
    rows = [
        {"gene": "G", "tissue": "Lung", "donor_id": "d1", "log2_tpm": 2.0},
        {"gene": "G", "tissue": "Lung", "donor_id": "d1", "log2_tpm": 4.0},
        {"gene": "G", "tissue": "Liver", "donor_id": "d2", "log2_tpm": 5.0},
    ]
    result = summarize_pan_tissue(rows, tpm_threshold=1.0)
    gene = result["genes"][0]
    assert gene["top_tissues"][0]["tissue"] == "Liver"
    assert gene["tissues"][1]["median"] == 3.0
    assert gene["tissues_with_majority_above_threshold"] == 2


def test_cell_type_expression_ranks_donor_level_detection():
    rows = [
        {
            "gene": "EPCAM",
            "cell_type": "epithelial",
            "cell_type_ontology_term_id": "CL:1",
            "donor_id": "A",
            "cells": 10,
            "fraction_detected": 0.8,
            "mean_log1p_raw_count": 1.2,
        },
        {
            "gene": "EPCAM",
            "cell_type": "epithelial",
            "cell_type_ontology_term_id": "CL:1",
            "donor_id": "B",
            "cells": 20,
            "fraction_detected": 0.6,
            "mean_log1p_raw_count": 1.0,
        },
        {
            "gene": "EPCAM",
            "cell_type": "immune",
            "cell_type_ontology_term_id": "CL:2",
            "donor_id": "A",
            "cells": 100,
            "fraction_detected": 0.1,
            "mean_log1p_raw_count": 0.2,
        },
    ]
    result = summarize_cell_type_expression(rows)
    epithelial = result["genes"][0]["top_cell_types"][0]
    assert epithelial["cell_type"] == "epithelial"
    assert epithelial["donors"] == 2
    assert epithelial["cells"] == 30
    assert epithelial["median_fraction_detected"] == 0.7
    assert epithelial["support"] == "limited"
    assert epithelial["leave_one_donor_out_detection_min"] == 0.6
    assert epithelial["leave_one_donor_out_detection_max"] == 0.8
    assert result["genes"][0]["leave_one_donor_out_top_agreement"] == 1.0
