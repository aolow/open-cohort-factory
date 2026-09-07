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

SAMPLES_QUERY = """(fn [dataset limit]
 (map :value (query {:select [:value] :from [:dataset]
  :join [:field [:= :dataset.id :dataset_id] :code [:= :field.id :field_id]]
  :limit limit :where [:and [:= :dataset.name dataset] [:= :field.name "sampleID"]]})))"""

PHENOTYPE_QUERY = """(fn [dataset samples fields]
 [nil (fetch [{:table dataset :columns fields :samples samples}])])"""

FIELD_CODES_QUERY = """(fn [dataset fields]
 (query {:select [:P.name [#sql/call
  [:group_concat :value :order :ordering :separator #sql/call [:chr 9]] :code]]
  :from [[{:select [:field.id :field.name] :from [:field]
   :join [{:table [[[:name :varchar fields]] :T]} [:= :T.name :field.name]]
   :where [:= :dataset_id {:select [:id] :from [:dataset]
    :where [:= :name dataset]}]} :P]]
  :left-join [:code [:= :P.id :field_id]] :group-by [:P.id]}))"""

PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"
PHENOTYPE_FIELDS = ["_study", "_primary_site", "_sample_type"]


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
                        "reference_context": (
                            sample.reference_context.value if sample.reference_context else None
                        ),
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

    def fetch_gtex_atlas(
        self, spec: ExpressionSpec
    ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
        """Fetch selected genes across every GTEx tissue represented in the Toil compendium."""
        samples = self._post(SAMPLES_QUERY, [PHENOTYPE_DATASET, None])
        _, encoded_columns = self._post(
            PHENOTYPE_QUERY, [PHENOTYPE_DATASET, samples, PHENOTYPE_FIELDS]
        )
        code_rows = self._post(FIELD_CODES_QUERY, [PHENOTYPE_DATASET, PHENOTYPE_FIELDS])
        codes = {row["name"]: row["code"].split("\t") for row in code_rows}
        decoded = {
            field: [_decode_code(value, codes[field]) for value in values]
            for field, values in zip(PHENOTYPE_FIELDS, encoded_columns, strict=True)
        }
        selected = [
            (sample, decoded["_primary_site"][index])
            for index, sample in enumerate(samples)
            if decoded["_study"][index] == "GTEX"
            and decoded["_sample_type"][index] == "Normal Tissue"
            and decoded["_primary_site"][index]
        ]
        identifiers = [sample for sample, _ in selected]
        tissue_by_sample = dict(selected)
        genes = list(dict.fromkeys(gene.strip().upper() for gene in spec.genes if gene.strip()))
        payload = self._post(GENE_QUERY, [spec.dataset, identifiers, genes])
        rows: list[dict[str, Any]] = []
        for result in payload:
            scores = result.get("scores") or []
            if len(scores) == 1 and isinstance(scores[0], list):
                scores = scores[0]
            for identifier, raw_value in zip(identifiers, scores, strict=False):
                value = _number_or_none(raw_value)
                if value is None:
                    continue
                tissue = tissue_by_sample[identifier]
                rows.append(
                    {
                        "sample_id": identifier,
                        "donor_id": "-".join(identifier.split("-")[:2]),
                        "cohort_role": "reference",
                        "cohort_name": f"GTEx pan-tissue: {tissue}",
                        "source": DataSource.GTEX.value,
                        "project_id": "GTEx Toil compendium",
                        "gene": str(result["gene"]).upper(),
                        "log2_tpm": _to_log2_tpm_plus_1(value),
                        "tissue": tissue,
                        "reference_context": "postmortem_reference",
                    }
                )
        return rows, {
            "phenotype_dataset": PHENOTYPE_DATASET,
            "eligible_gtex_samples": len(identifiers),
            "returned_gtex_samples": len({row["sample_id"] for row in rows}),
            "tissue_count": len({row["tissue"] for row in rows}),
        }

    def _post(self, query_function: str, arguments: list[Any]) -> Any:
        encoded = ["nil" if argument is None else json.dumps(argument) for argument in arguments]
        response = self._client.post(
            f"{self.host}/data/",
            content=f"({query_function} {' '.join(encoded)})",
            headers={"Content-Type": "text/plain"},
        )
        response.raise_for_status()
        return response.json()


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


def _decode_code(value: Any, codes: list[str]) -> str | None:
    number = _number_or_none(value)
    if number is None:
        return None
    index = int(number)
    return codes[index] if 0 <= index < len(codes) else None
