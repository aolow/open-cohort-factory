from open_cohort_factory.expression import summarize_expression


def test_expression_reports_unmatched_and_donor_matched_results():
    rows = [
        {"gene": "EPCAM", "cohort_role": "disease", "donor_id": "d1", "log2_tpm": 8.0},
        {"gene": "EPCAM", "cohort_role": "disease", "donor_id": "d2", "log2_tpm": 4.0},
        {"gene": "EPCAM", "cohort_role": "reference", "donor_id": "r1", "log2_tpm": 2.0},
        {"gene": "EPCAM", "cohort_role": "reference", "donor_id": "r2", "log2_tpm": 0.0},
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
