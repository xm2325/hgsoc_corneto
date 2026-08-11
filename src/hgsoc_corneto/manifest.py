"""Construction of the normalized OCM master manifest."""

from __future__ import annotations

from collections import defaultdict
from typing import Any


def _passage_number(value: str | None) -> int:
    if not value:
        return -1
    return int(value.removeprefix("P"))


def build_master_manifest(
    rna_runs: list[dict[str, Any]],
    tighe_rows: list[dict[str, Any]],
    abcb1_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Join source-derived metadata and assign label-blind eligibility flags."""

    tighe_by_ocm = {row["canonical_ocm_id"]: row for row in tighe_rows}
    abcb1_by_ocm = {row["canonical_ocm_id"]: row for row in abcb1_rows}
    if set(tighe_by_ocm) != set(abcb1_by_ocm):
        only_s1 = sorted(set(tighe_by_ocm) - set(abcb1_by_ocm))
        only_s2 = sorted(set(abcb1_by_ocm) - set(tighe_by_ocm))
        raise ValueError(f"Tighe S1/S2 OCM mismatch: S1-only={only_s1}; S2-only={only_s2}")

    groups: dict[tuple[str | None, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rna_runs:
        groups[(row["canonical_ocm_id"], row["sample_class"])].append(row)

    representative_runs: set[str] = set()
    for (ocm_id, _sample_class), rows in groups.items():
        if ocm_id is None:
            continue
        chosen = sorted(
            rows,
            key=lambda row: (
                _passage_number(row.get("passage")),
                row["study_accession"],
                row["run_accession"],
            ),
        )[0]
        representative_runs.add(chosen["run_accession"])

    output: list[dict[str, Any]] = []
    for run in rna_runs:
        ocm_id = run["canonical_ocm_id"]
        tighe = tighe_by_ocm.get(ocm_id, {})
        abcb1 = abcb1_by_ocm.get(ocm_id, {})
        representative = run["run_accession"] in representative_runs
        hgsoc_tumour = (
            run["sample_class"] == "tumour"
            and tighe.get("histotype_group") == "HGSOC"
            and representative
        )
        record = {
            **run,
            "in_tighe_83_ocm_screen": bool(tighe),
            "histotype_reported": tighe.get("histotype_reported"),
            "histotype_group": tighe.get("histotype_group"),
            "figo_stage": tighe.get("figo_stage"),
            "chemo_naive_at_biopsy": tighe.get("chemo_naive_at_biopsy"),
            "biopsy_type": tighe.get("biopsy_type"),
            "tighe_table_row": tighe.get("table_row"),
            "abcb1_normalized_read_count": abcb1.get("abcb1_normalized_read_count"),
            "abcb1_missing_reason": abcb1.get("abcb1_missing_reason"),
            "is_representative_rna_library": representative,
            "hgsoc_tumour_eligible": hgsoc_tumour,
            "primary_cohort_eligible": hgsoc_tumour and bool(tighe),
            "exact_paclitaxel_auc_available": False,
            "exact_paclitaxel_gi50_available": False,
            "exact_cumulative_paclitaxel_exposure_available": False,
            "paclitaxel_auc": None,
            "paclitaxel_gi50_nm": None,
            "cumulative_paclitaxel_exposure_mg": None,
        }
        output.append(record)
    return sorted(output, key=lambda row: (row["study_accession"], row["source_name"]))
