"""Turn a cohort build summary into a concise analysis decision memo."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def decision_memo(summary: dict[str, Any]) -> str:
    lines = ["# Cohort comparability decision", ""]
    for comparison in summary.get("comparability", []):
        female = comparison["female_proportion"]
        lines.extend(
            [
                f"## {comparison['reference_panel']}",
                "",
                f"**Decision:** `{comparison['decision']}`",
                "",
                f"- Disease donors: {comparison['disease_donors']}",
                f"- Reference donors: {comparison['reference_donors']}",
                "- Disease age brackets: "
                f"{json.dumps(comparison['disease_age_brackets'], sort_keys=True)}",
                "- Reference age brackets: "
                f"{json.dumps(comparison['reference_age_brackets'], sort_keys=True)}",
                "- Uncovered disease ages: "
                f"{', '.join(comparison['uncovered_disease_age_brackets']) or 'none'}",
                f"- Female proportion: disease {female['disease']:.1%}, "
                f"reference {female['reference']:.1%}",
                "",
                "### Required analysis actions",
                "",
            ]
        )
        lines.extend(
            f"{index}. {action}"
            for index, action in enumerate(comparison["recommended_actions"], 1)
        )
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="Path to a generated summary.json")
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    print(decision_memo(summary))


if __name__ == "__main__":
    main()
