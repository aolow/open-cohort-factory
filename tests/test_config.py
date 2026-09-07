from pathlib import Path

import pytest
from pydantic import ValidationError

from open_cohort_factory.config import load_spec
from open_cohort_factory.models import AgeRange, ReferenceContext


def test_example_spec_is_valid() -> None:
    spec = load_spec(Path("examples/luad_older_adults.yaml"))
    assert spec.disease_cohort.projects == ["TCGA-LUAD"]
    assert spec.reference_panels[1].context == ReferenceContext.POSTMORTEM_REFERENCE
    assert spec.expression is not None and spec.expression.gtex_pan_tissue


def test_age_range_must_be_ordered() -> None:
    with pytest.raises(ValidationError):
        AgeRange(minimum=70, maximum=60)
