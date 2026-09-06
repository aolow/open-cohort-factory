"""Small, provenance-aware client for the UCSC Xena Toil expression compendium."""

from __future__ import annotations

import json
import math
import re
from datetime import UTC, datetime
from typing import Any

import httpx

from ..models import DataSource, ExpressionSpec, ProvenanceRecord, SampleRecord

XENA_HOST = "https://toil.xenahubs.net"
XENA_CITATION_URL = (
    "https://xenabrowser.net/datapages/"
    "?cohort=tcga+target+gtex&hub=https%3A%2F%2Ftoil.xenahubs.net%3A443"
)

# This is the official xenaPython datasetGeneProbeAvg query. Keeping the query here
# avoids adding an unmaintained runtime dependency while using Xena's public API.
GENE_QUERY = """(fn [dataset samples genes]
 (let [probemap (:probemap (car (query {:select [:probemap] :from [:dataset]
   :where [:= :name dataset]})))
  get-probes (fn [gene] (xena-query {:select ["name" "position"] :from [probemap]
   :where [:in :any "genes" [gene]]}))
  avg (fn [scores] (mean scores 0))
  scores-for-gene (fn [gene]
   (let [probes (get-probes gene) probe-names (probes "name")
    scores (fetch [{:table dataset :samples samples :columns probe-names}])]
    {:gene gene :position (probes "position")
     :scores (if (car probe-names) (avg scores) [[]])}))]
  (map scores-for-gene genes)))"""


def xena_sample_id(sample: SampleRecord) -> str | None:
    """Return the public Toil sample identifier without guessing from UUIDs."""
    if sample.source == DataSource.GDC:
        identifier = sample.sample_submitter_id
        if identifier and re.fullmatch(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}[A-Z].*", identifier):
            # Toil indexes TCGA at sample level; GDC exposes the longer vial/portion barcode.
            return identifier[:15]
        return identifier
    if sample.source == DataSource.GTEX:
        return sample.sample_id
    return None


class XenaClient:
    """Fetch small gene slices from the uniformly processed Toil matrix."""

    def __init__(self, client: httpx.Client | None = None, host: str = XENA_HOST) -> None:
        self._client = client or httpx.Client(timeout=120, follow_redirects=True)
        self.host = host.rstrip("/")

    def fetch(
        self, samples: list[SampleRecord], spec: ExpressionSpec
    ) -> tuple[list[dict[str, Any]], ProvenanceRecord]:
        eligible = [(sample, xena_sample_id(sample)) for sample in samples]
        eligible = [(sample, identifier) for sample, identifier in eligible if identifier]
        identifiers = list(dict.fromkeys(identifier for _, identifier in eligible))
        genes = list(dict.fromkeys(gene.strip().upper() for gene in spec.genes if gene.strip()))
        if not identifiers:
            raise ValueError("No cohort samples have public Xena identifiers.")

        arguments = " ".join(
            [json.dumps(spec.dataset), json.dumps(identifiers), json.dumps(genes)]
        )
        query = f"({GENE_QUERY} {arguments})"
        response = self._client.post(
            f"{self.host}/data/", content=query, headers={"Content-Type": "text/plain"}
        )
        response.raise_for_status()
        payload = response.json()
        by_id = {identifier: sample for sample, identifier in eligible}
        rows: list[dict[str, Any]] = []
        for result in payload:
            gene = str(result["gene"]).upper()
            scores = result.get("scores") or []
            # Xena's mean() result is wrapped once (one vector per requested gene).
            if len(scores) == 1 and isinstance(scores[0], list):
                scores = scores[0]
            for identifier, raw_value in zip(identifiers, scores, strict=False):
                value = _number_or_none(raw_value)
                if value is None:
                    continue
                sample = by_id[identifier]
                rows.append(
                    {
                        "sample_id": identifier,
                        "donor_id": sample.case_id,
                        "cohort_role": sample.cohort_role,
                        "cohort_name": sample.cohort_name,
                        "source": sample.source.value,
                        "project_id": sample.project_id,
                        "gene": gene,
                        "log2_tpm": _to_log2_tpm_plus_1(value),
                    }
                )

        provenance = ProvenanceRecord(
            source=DataSource.XENA_TOIL,
            endpoint=f"{self.host}/data/",
            retrieved_at=datetime.now(UTC).isoformat(),
            query={
                "dataset": spec.dataset,
                "genes": genes,
                "requested_sample_count": len(identifiers),
                "returned_sample_count": len({row["sample_id"] for row in rows}),
                "returned_samples_by_source": {
                    source: len({row["sample_id"] for row in rows if row["source"] == source})
                    for source in sorted({row["source"] for row in rows})
                },
                "transform": spec.transform,
                "input_scale": "log2(TPM + 0.001)",
            },
            source_release="hg38 / Gencode v23 / Toil RSEM recompute",
            citation_url=XENA_CITATION_URL,
        )
        return rows, provenance


def _number_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _to_log2_tpm_plus_1(value: float) -> float:
    """Convert Xena log2(TPM + 0.001) values to interpretable log2(TPM + 1)."""
    tpm = max(0.0, 2**value - 0.001)
    return math.log2(tpm + 1.0)
