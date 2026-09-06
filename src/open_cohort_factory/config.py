"""Configuration loading."""

from pathlib import Path

import yaml

from .models import ProjectSpec


def load_spec(path: Path) -> ProjectSpec:
    """Load and validate a cohort specification from YAML."""
    with path.open(encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    return ProjectSpec.model_validate(payload)

