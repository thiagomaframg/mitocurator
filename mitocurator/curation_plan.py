from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple

from .utils import ensure_dir


def _safe_read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, "r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _contains_any(text: str, terms: List[str]) -> bool:
    text = (text or "").upper()
    return any(term.upper() in text for term in terms)


def _infer_downstream_actions(row: Dict[str, str]) -> Tuple[List[str], str]:
    gene = row.get("gene", ".")
    priority = row.get("annotation_priority", ".")
    issues = row.get("annotation_issues", "")
    recommendations = row.get("annotation_recommendations", "")
    interpretation = row.get("annotation_interpretation", "")

    combined = " ".join([gene, priority, issues, recommendations, interpretation]).upper()
    actions: List[str] = []
    reasons: List[str] = []

    # CDS with internal stops: first downstream evidence should usually be read support.
    if _contains_any(combined, ["INTERNAL_STOP", "CHECK_READ_SUPPORT", "STOP"]):
        actions.append("read_support")
        reasons.append("CDS has internal-stop or stop-codon evidence; validate with read support")

    # Missing/partial CDS with candidate evidence: targeted extraction/consensus path.
    if _contains_any(combined, ["MISSING_GENE", "MISSING", "PARTIAL_REF_MATCH", "STRONG_CANDIDATE_REGION", "REVIEW_PARTIAL_CANDIDATE"]):
        if not gene.upper().startswith("TRNA"):
            actions.extend([
                "targeted_extraction",
                "reconstruction_pools",
                "targeted_consensus",
            ])
            reasons.append("Expected CDS is missing/partial or has candidate evidence; use targeted consensus workflow")

    # Candidate assembly is useful after targeted consensus for difficult missing/partial CDS cases.
    if _contains_any(combined, ["MISSING_GENE", "STRONG_CANDIDATE_REGION", "PARTIAL_REF_MATCH"]) and not gene.upper().startswith("TRNA"):
        actions.append("candidate_assembly")
        reasons.append("Missing/partial CDS may require local assembly validation")

    # Missing tRNAs should not trigger heavy read-backed CDS workflows by default.
    if gene.upper().startswith("TRNA") and _contains_any(combined, ["MISSING", "REVIEW_MISSING_EXPECTED_GENE"]):
        actions.append("manual_trna_review")
        reasons.append("Expected tRNA is missing; review tRNA caller output and naming/anticodon conventions")

    if not actions:
        actions.append("manual_review")
        reasons.append("No deterministic downstream action inferred; manual review recommended")

    # Preserve order and remove duplicates.
    unique_actions = []
    for action in actions:
        if action not in unique_actions:
            unique_actions.append(action)

    return unique_actions, "; ".join(reasons)


def build_downstream_curation_plan(root: Path) -> Tuple[Path, Dict[str, bool]]:
    """Build downstream curation plan from annotation review targets.

    The plan is diagnostic/planning output. It does not modify GenBank files.
    It indicates which downstream stages should be considered or auto-run by
    the pipeline based on annotation-level evidence.
    """
    assessment_dir = root / "06_annotation_assessment"
    review_targets = assessment_dir / "annotation_review_targets.tsv"
    plan_tsv = assessment_dir / "downstream_curation_plan.tsv"

    rows = _safe_read_tsv(review_targets)
    plan_rows: List[Dict[str, str]] = []

    run_flags = {
        "read_support": False,
        "targeted_extraction": False,
        "reconstruction_pools": False,
        "targeted_consensus": False,
        "candidate_assembly": False,
    }

    for row in rows:
        actions, reason = _infer_downstream_actions(row)
        for action in actions:
            if action in run_flags:
                run_flags[action] = True

        plan_rows.append({
            "gene": row.get("gene", "."),
            "annotation_priority": row.get("annotation_priority", "."),
            "annotation_issues": row.get("annotation_issues", "."),
            "annotation_recommendations": row.get("annotation_recommendations", "."),
            "recommended_downstream_actions": ";".join(actions),
            "plan_reason": reason,
        })

    _write_tsv(
        plan_tsv,
        plan_rows,
        [
            "gene",
            "annotation_priority",
            "annotation_issues",
            "annotation_recommendations",
            "recommended_downstream_actions",
            "plan_reason",
        ],
    )

    return plan_tsv, run_flags
