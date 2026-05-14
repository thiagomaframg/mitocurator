from __future__ import annotations

import csv
import html
from pathlib import Path
from typing import Dict, List

from .utils import ensure_dir


ANNOTATION_EVIDENCE_PATHS = [
    "05_refinement/refined.gb",
    "05_refinement/expected_gene_set.tsv",
    "05_refinement/added_features.tsv",
    "05_refinement/missing_gene_candidates.tsv",
    "05_refinement/cds_refinement_candidates.tsv",
    "05_refinement/reference_similarity_candidates.tsv",
    "05_refinement/problematic_cds_reference_check.tsv",
    "05_refinement/problematic_cds_stop_context.tsv",
    "05_refinement/problematic_cds_reference_alignment.tsv",
    "05_refinement/curation_recommendations.tsv",
    "05_refinement/curation_recommendations.md",
    "07_gene_qc/gene_qc.tsv",
    "07_gene_qc/problematic_features.tsv",
    "07_gene_qc/intergenic_regions.tsv",
    "07_gene_qc/sequence_summary.tsv",
    "07_gene_qc/diagnostic_summary.md",
]


PRIORITY_RANK = {
    "HIGH": 3,
    "MEDIUM": 2,
    "LOW": 1,
    "REVIEW": 1,
    "UNKNOWN": 0,
    ".": 0,
    "": 0,
}


def _safe_read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except Exception:
        return []


def _first_existing_column(fieldnames: List[str], candidates: List[str]) -> str | None:
    for c in candidates:
        if c in fieldnames:
            return c
    return None


def _collect_evidence_status(root: Path) -> List[Dict[str, str]]:
    rows = []
    for rel in ANNOTATION_EVIDENCE_PATHS:
        p = root / rel
        rows.append({
            "file": rel,
            "status": "available" if p.exists() else "not available",
        })
    return rows


def _normalize_priority(value: str | None, recommendation: str | None = None) -> str:
    raw = str(value or "").strip().upper()
    if raw in {"HIGH", "MEDIUM", "LOW"}:
        return raw

    rec = str(recommendation or raw).strip().upper()

    high_terms = [
        "CHECK_READ_SUPPORT_FOR_LOCAL_ERROR",
        "CHECK_FRAMESHIFT_OR_LOCAL_INDEL",
        "REVIEW_FOR_MANUAL_ANNOTATION",
        "CHECK_INTERNAL_STOP",
        "MISSING_GENE",
        "INTERNAL_STOP",
    ]
    medium_terms = [
        "REVIEW_PARTIAL_CANDIDATE",
        "SEARCH_EXTENDED_REGION_OR_REANNOTATE_BOUNDARIES",
        "MANUAL_REVIEW",
        "PARTIAL_REF_MATCH",
    ]
    low_terms = [
        "NO_CANDIDATE_FOUND",
        "POSSIBLE_MISANNOTATION",
        "WEAK_REF_MATCH",
    ]

    if any(term in rec for term in high_terms):
        return "HIGH"
    if any(term in rec for term in medium_terms):
        return "MEDIUM"
    if any(term in rec for term in low_terms):
        return "LOW"

    return "REVIEW"


def _append_row(
    rows: List[Dict[str, str]],
    gene: str,
    issue_type: str,
    status: str,
    priority: str | None,
    recommendation: str,
    evidence_source: str,
    evidence_summary: str,
):
    recommendation = str(recommendation or ".").strip() or "."
    rows.append({
        "gene": str(gene or ".").strip() or ".",
        "issue_type": str(issue_type or ".").strip() or ".",
        "status": str(status or ".").strip() or ".",
        "priority": _normalize_priority(priority, recommendation),
        "recommendation": recommendation,
        "evidence_source": evidence_source,
        "evidence_summary": str(evidence_summary or ".").strip() or ".",
    })


def _collect_annotation_rows(root: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []

    # 1. Expected gene set: present/missing genes
    source = "05_refinement/expected_gene_set.tsv"
    data = _safe_read_tsv(root / source)
    if data:
        cols = list(data[0].keys())
        gene_col = _first_existing_column(cols, ["gene", "target"])
        type_col = _first_existing_column(cols, ["type", "target_type"])
        status_col = _first_existing_column(cols, ["status"])
        count_col = _first_existing_column(cols, ["count", "n", "observed_count"])

        for r in data:
            gene = r.get(gene_col, ".") if gene_col else "."
            target_type = r.get(type_col, ".") if type_col else "."
            status = r.get(status_col, ".") if status_col else "."
            count = r.get(count_col, ".") if count_col else "."
            issue = "expected_gene_status"
            priority = "HIGH" if status.upper() == "MISSING" and target_type.upper() == "CDS" else "LOW"
            recommendation = "REVIEW_MISSING_EXPECTED_GENE" if status.upper() == "MISSING" else "NO_ACTION_EXPECTED_GENE_PRESENT"
            _append_row(
                rows,
                gene=gene,
                issue_type=issue,
                status=status,
                priority=priority,
                recommendation=recommendation,
                evidence_source=source,
                evidence_summary=f"type={target_type}; count={count}",
            )

    # 2. Problematic features from gene QC
    source = "07_gene_qc/problematic_features.tsv"
    data = _safe_read_tsv(root / source)
    if data:
        cols = list(data[0].keys())
        gene_col = _first_existing_column(cols, ["gene_normalized", "gene", "feature", "feature_name"])
        status_col = _first_existing_column(cols, ["decision_hint", "status", "comment"])
        feature_col = _first_existing_column(cols, ["type", "feature_type", "kind"])

        for r in data:
            hint = r.get(status_col, "CHECK_FEATURE") if status_col else "CHECK_FEATURE"
            feature_type = r.get(feature_col, ".") if feature_col else "."
            _append_row(
                rows,
                gene=r.get(gene_col, ".") if gene_col else ".",
                issue_type="problematic_feature",
                status="PROBLEMATIC",
                priority="HIGH",
                recommendation=hint,
                evidence_source=source,
                evidence_summary=f"feature_type={feature_type}; hint={hint}",
            )

    # 3. Curation recommendations from refinement
    source = "05_refinement/curation_recommendations.tsv"
    data = _safe_read_tsv(root / source)
    if data:
        cols = list(data[0].keys())
        gene_col = _first_existing_column(cols, ["target", "gene", "gene_normalized", "target_gene"])
        issue_col = _first_existing_column(cols, ["issue_type", "problem", "issue", "category"])
        status_col = _first_existing_column(cols, ["status"])
        priority_col = _first_existing_column(cols, ["priority"])
        rec_col = _first_existing_column(cols, ["recommendation"])
        evidence_col = _first_existing_column(cols, ["evidence_summary", "comment"])

        for r in data:
            recommendation = r.get(rec_col, ".") if rec_col else "."
            _append_row(
                rows,
                gene=r.get(gene_col, ".") if gene_col else ".",
                issue_type=r.get(issue_col, "refinement_recommendation") if issue_col else "refinement_recommendation",
                status=r.get(status_col, ".") if status_col else ".",
                priority=r.get(priority_col, None) if priority_col else None,
                recommendation=recommendation,
                evidence_source=source,
                evidence_summary=r.get(evidence_col, recommendation) if evidence_col else recommendation,
            )

    # 4. Missing gene candidates
    source = "05_refinement/missing_gene_candidates.tsv"
    data = _safe_read_tsv(root / source)
    if data:
        cols = list(data[0].keys())
        gene_col = _first_existing_column(cols, ["gene"])
        hint_col = _first_existing_column(cols, ["decision_hint", "comment"])
        cand_col = _first_existing_column(cols, ["candidate_id"])
        length_col = _first_existing_column(cols, ["length_aa", "candidate_aa_length"])

        for r in data:
            hint = r.get(hint_col, ".") if hint_col else "."
            _append_row(
                rows,
                gene=r.get(gene_col, ".") if gene_col else ".",
                issue_type="missing_gene_candidate",
                status="CANDIDATE",
                priority=_normalize_priority(hint),
                recommendation=hint,
                evidence_source=source,
                evidence_summary=f"candidate={r.get(cand_col, '.') if cand_col else '.'}; length_aa={r.get(length_col, '.') if length_col else '.'}",
            )

    # 5. CDS refinement candidates
    source = "05_refinement/cds_refinement_candidates.tsv"
    data = _safe_read_tsv(root / source)
    if data:
        cols = list(data[0].keys())
        gene_col = _first_existing_column(cols, ["gene"])
        hint_col = _first_existing_column(cols, ["decision_hint", "comment"])
        old_stop_col = _first_existing_column(cols, ["old_internal_stop_count"])
        new_stop_col = _first_existing_column(cols, ["candidate_internal_stop_count"])

        for r in data:
            hint = r.get(hint_col, ".") if hint_col else "."
            _append_row(
                rows,
                gene=r.get(gene_col, ".") if gene_col else ".",
                issue_type="cds_refinement_candidate",
                status="CANDIDATE",
                priority=_normalize_priority(hint),
                recommendation=hint,
                evidence_source=source,
                evidence_summary=(
                    f"old_internal_stop_count={r.get(old_stop_col, '.') if old_stop_col else '.'}; "
                    f"candidate_internal_stop_count={r.get(new_stop_col, '.') if new_stop_col else '.'}"
                ),
            )

    # 6. Reference similarity candidates
    source = "05_refinement/reference_similarity_candidates.tsv"
    data = _safe_read_tsv(root / source)
    if data:
        cols = list(data[0].keys())
        gene_col = _first_existing_column(cols, ["gene"])
        hint_col = _first_existing_column(cols, ["decision_hint"])
        pid_col = _first_existing_column(cols, ["percent_identity"])
        cov_col = _first_existing_column(cols, ["aligned_coverage_reference"])
        ratio_col = _first_existing_column(cols, ["length_ratio"])

        for r in data:
            hint = r.get(hint_col, ".") if hint_col else "."
            _append_row(
                rows,
                gene=r.get(gene_col, ".") if gene_col else ".",
                issue_type="reference_similarity_candidate",
                status="CANDIDATE",
                priority=_normalize_priority(hint),
                recommendation=hint,
                evidence_source=source,
                evidence_summary=(
                    f"pid={r.get(pid_col, '.') if pid_col else '.'}; "
                    f"cov_ref={r.get(cov_col, '.') if cov_col else '.'}; "
                    f"length_ratio={r.get(ratio_col, '.') if ratio_col else '.'}"
                ),
            )

    return rows


def _collapse_annotation_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    grouped: Dict[tuple, Dict[str, object]] = {}

    for r in rows:
        key = (
            r["gene"],
            r["issue_type"],
            r["status"],
            r["priority"],
            r["recommendation"],
            r["evidence_source"],
        )
        if key not in grouped:
            grouped[key] = {
                "gene": r["gene"],
                "issue_type": r["issue_type"],
                "status": r["status"],
                "priority": r["priority"],
                "recommendation": r["recommendation"],
                "evidence_source": r["evidence_source"],
                "n_records": 0,
                "evidence_examples": [],
            }

        grouped[key]["n_records"] = int(grouped[key]["n_records"]) + 1
        if r["evidence_summary"] not in grouped[key]["evidence_examples"]:
            grouped[key]["evidence_examples"].append(r["evidence_summary"])

    collapsed = []
    for item in grouped.values():
        examples = item["evidence_examples"]
        evidence_examples = "; ".join(examples[:4])
        if len(examples) > 4:
            evidence_examples += f"; ... (+{len(examples) - 4} more)"

        collapsed.append({
            "gene": str(item["gene"]),
            "issue_type": str(item["issue_type"]),
            "status": str(item["status"]),
            "priority": str(item["priority"]),
            "recommendation": str(item["recommendation"]),
            "evidence_source": str(item["evidence_source"]),
            "n_records": str(item["n_records"]),
            "evidence_examples": evidence_examples or ".",
        })

    collapsed.sort(
        key=lambda r: (
            -PRIORITY_RANK.get(r["priority"], 0),
            r["gene"],
            r["issue_type"],
            r["recommendation"],
        )
    )
    return collapsed


def _build_gene_annotation_summary(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    # The gene-level summary is intended to guide manual review.
    # Therefore, annotation rows explicitly classified as "no action" are kept
    # in the detailed table but excluded from the decision-oriented summary.
    review_rows = [
        r for r in rows
        if r.get("recommendation") != "NO_ACTION_EXPECTED_GENE_PRESENT"
    ]

    by_gene: Dict[str, List[Dict[str, str]]] = {}
    for r in review_rows:
        by_gene.setdefault(r["gene"], []).append(r)

    output = []
    for gene, gene_rows in sorted(by_gene.items()):
        final_priority = max(
            [r["priority"] for r in gene_rows],
            key=lambda p: PRIORITY_RANK.get(p, 0),
        )

        issues = sorted({r["issue_type"] for r in gene_rows})
        recommendations = sorted({r["recommendation"] for r in gene_rows})
        evidence_sources = sorted({r["evidence_source"] for r in gene_rows})

        has_missing = any("MISSING" in r["issue_type"].upper() or "MISSING" in r["recommendation"].upper() for r in gene_rows)
        has_stop = any("STOP" in r["issue_type"].upper() or "STOP" in r["recommendation"].upper() or "INTERNAL_STOP" in r["issue_type"].upper() for r in gene_rows)
        has_candidate = any("candidate" in r["issue_type"].lower() for r in gene_rows)

        is_trna = gene.lower().startswith("trna")

        if has_missing and is_trna:
            interpretation = "Expected tRNA is missing; review tRNA caller output, anticodon naming and possible synonym conventions."
        elif final_priority == "HIGH" and has_missing:
            interpretation = "Expected CDS is missing or only represented by candidate evidence; manual annotation review is required before downstream read-backed curation."
        elif final_priority == "HIGH" and has_stop:
            interpretation = "Annotated CDS has internal-stop evidence; inspect boundaries, frame and reference similarity before downstream read-backed curation."
        elif final_priority == "HIGH":
            interpretation = "High-priority annotation issue detected; manual review is required before downstream curation."
        elif final_priority == "MEDIUM" and has_candidate:
            interpretation = "Candidate annotation evidence exists but is partial or ambiguous; review before accepting."
        elif final_priority == "MEDIUM":
            interpretation = "Medium-priority annotation review recommended."
        elif final_priority == "LOW":
            interpretation = "Low-priority annotation signal; review if biologically relevant or if expected features are missing."
        else:
            interpretation = "Annotation evidence available, but priority could not be fully classified."

        output.append({
            "gene": gene,
            "annotation_priority": final_priority,
            "n_annotation_evidence_sources": str(len(evidence_sources)),
            "annotation_evidence_sources": ";".join(evidence_sources),
            "annotation_issues": ";".join(issues),
            "annotation_recommendations": ";".join(recommendations),
            "annotation_interpretation": interpretation,
        })

    output.sort(key=lambda r: (-PRIORITY_RANK.get(r["annotation_priority"], 0), r["gene"]))
    return output


def _build_gene_inventory(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Build a complete annotation inventory from annotation-level evidence.

    This table is intended to include both apparently OK genes and genes that
    require review. It is complementary to annotation_review_targets.tsv, which
    is decision-oriented and excludes explicit no-action rows.
    """
    by_gene: Dict[str, List[Dict[str, str]]] = {}
    for r in rows:
        by_gene.setdefault(r["gene"], []).append(r)

    inventory = []
    for gene, gene_rows in sorted(by_gene.items()):
        priorities = [r["priority"] for r in gene_rows]
        annotation_priority = max(
            priorities,
            key=lambda p: PRIORITY_RANK.get(p, 0),
        ) if priorities else "UNKNOWN"

        issues = sorted({r["issue_type"] for r in gene_rows})
        recommendations = sorted({r["recommendation"] for r in gene_rows})
        evidence_sources = sorted({r["evidence_source"] for r in gene_rows})

        expected_rows = [r for r in gene_rows if r["issue_type"] == "expected_gene_status"]
        expected_status = "."
        gene_type = "."
        if expected_rows:
            # expected_gene_status rows store type/count in evidence_examples, e.g. type=CDS; count=1
            ex = expected_rows[0].get("evidence_examples", "")
            expected_status = expected_rows[0].get("status", ".")
            for part in ex.split(";"):
                part = part.strip()
                if part.startswith("type="):
                    gene_type = part.replace("type=", "", 1).strip() or "."

        has_no_action_only = all(
            r["recommendation"] == "NO_ACTION_EXPECTED_GENE_PRESENT"
            for r in gene_rows
        )

        has_missing = any(
            "MISSING" in r["status"].upper()
            or "MISSING" in r["issue_type"].upper()
            or "MISSING" in r["recommendation"].upper()
            for r in gene_rows
        )
        has_problem = any(
            r["status"].upper() == "PROBLEMATIC"
            or "STOP" in r["recommendation"].upper()
            or "INTERNAL_STOP" in r["issue_type"].upper()
            or "problematic" in r["issue_type"].lower()
            for r in gene_rows
        )
        has_candidate = any("candidate" in r["issue_type"].lower() for r in gene_rows)

        if has_no_action_only:
            annotation_status = "OK"
            recommended_action = "KEEP"
            interpretation = "Expected feature is present and no annotation-level issue was detected."
        elif has_missing:
            annotation_status = "MISSING_OR_INCOMPLETE"
            recommended_action = "REVIEW_MISSING_OR_CANDIDATE_FEATURE"
            interpretation = "Expected feature is missing or represented only by candidate evidence."
        elif has_problem:
            annotation_status = "PROBLEMATIC"
            recommended_action = "REVIEW_ANNOTATION"
            interpretation = "Annotation-level problem detected; inspect coordinates, frame and reference evidence."
        elif has_candidate:
            annotation_status = "CANDIDATE_REVIEW"
            recommended_action = "REVIEW_CANDIDATE"
            interpretation = "Candidate evidence exists; review before accepting the annotation."
        else:
            annotation_status = "REVIEW"
            recommended_action = "MANUAL_REVIEW"
            interpretation = "Annotation evidence exists but could not be classified as OK or problematic."

        inventory.append({
            "gene": gene,
            "type": gene_type,
            "expected_status": expected_status,
            "annotation_status": annotation_status,
            "annotation_priority": annotation_priority,
            "recommended_action": recommended_action,
            "annotation_issues": ";".join(issues),
            "annotation_recommendations": ";".join(recommendations),
            "annotation_evidence_sources": ";".join(evidence_sources),
            "annotation_interpretation": interpretation,
        })

    inventory.sort(
        key=lambda r: (
            0 if r["annotation_status"] != "OK" else 1,
            -PRIORITY_RANK.get(r["annotation_priority"], 0),
            r["type"],
            r["gene"],
        )
    )
    return inventory


def _html_escape(value: object) -> str:
    return html.escape(str(value if value is not None else "."))


def _status_class(value: str) -> str:
    v = str(value or "").upper()
    if v in {"OK", "KEEP", "PRESENT"}:
        return "ok"
    if "HIGH" in v or "PROBLEM" in v or "MISSING" in v:
        return "problem"
    if "MEDIUM" in v or "REVIEW" in v or "CANDIDATE" in v:
        return "warning"
    if "LOW" in v:
        return "low"
    return "neutral"


def _count_by(rows: List[Dict[str, str]], key: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for row in rows:
        value = str(row.get(key, ".") or ".")
        counts[value] = counts.get(value, 0) + 1
    return counts


def _html_escape(value: object) -> str:
    return html.escape(str(value if value is not None else "."))


def _status_class(value: str) -> str:
    v = str(value or "").upper()
    if v in {"OK", "KEEP", "PRESENT"}:
        return "ok"
    if "HIGH" in v or "PROBLEM" in v or "MISSING" in v:
        return "problem"
    if "MEDIUM" in v or "REVIEW" in v or "CANDIDATE" in v:
        return "warning"
    if "LOW" in v:
        return "low"
    return "neutral"


def _count_by(rows: List[Dict[str, str]], key: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for row in rows:
        value = str(row.get(key, ".") or ".")
        counts[value] = counts.get(value, 0) + 1
    return counts


def _write_html_report(
    html_path: Path,
    evidence: List[Dict[str, str]],
    inventory_rows: List[Dict[str, str]],
    review_rows: List[Dict[str, str]],
    evidence_rows: List[Dict[str, str]],
    md_path: Path,
    annotation_summary_tsv: Path,
    gene_inventory_tsv: Path,
    gene_summary_tsv: Path,
):
    total_genes = len(inventory_rows)
    ok_genes = sum(1 for r in inventory_rows if r.get("annotation_status") == "OK")
    problem_genes = sum(1 for r in inventory_rows if r.get("annotation_status") == "PROBLEMATIC")
    missing_genes = sum(1 for r in inventory_rows if r.get("annotation_status") == "MISSING_OR_INCOMPLETE")
    review_targets = len(review_rows)
    high_targets = sum(1 for r in review_rows if r.get("annotation_priority") == "HIGH")

    available_files = [r["file"] for r in evidence if r.get("status") == "available"]
    unavailable_files = [r["file"] for r in evidence if r.get("status") != "available"]

    by_type = _count_by(inventory_rows, "type")
    by_status = _count_by(inventory_rows, "annotation_status")

    css = """
:root {
  --bg: #f6f8fb;
  --panel: #ffffff;
  --text: #1f2937;
  --muted: #6b7280;
  --border: #e5e7eb;
  --ok: #047857;
  --ok-bg: #d1fae5;
  --warn: #b45309;
  --warn-bg: #fef3c7;
  --problem: #b91c1c;
  --problem-bg: #fee2e2;
  --low: #475569;
  --low-bg: #e2e8f0;
  --accent: #1d4ed8;
  --accent-bg: #dbeafe;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
  line-height: 1.45;
}
.container {
  max-width: 1180px;
  margin: 0 auto;
  padding: 32px 20px 56px;
}
.header {
  background: linear-gradient(135deg, #1e3a8a, #2563eb);
  color: white;
  border-radius: 20px;
  padding: 28px 32px;
  box-shadow: 0 12px 30px rgba(37, 99, 235, .18);
}
.header h1 { margin: 0 0 8px; font-size: 30px; }
.header p { margin: 0; opacity: .92; }

.grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 14px;
  margin: 20px 0;
}
.card {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 18px;
  padding: 18px;
  box-shadow: 0 8px 20px rgba(15, 23, 42, .05);
}
.card .label {
  color: var(--muted);
  font-size: 13px;
  margin-bottom: 8px;
}
.card .value {
  font-size: 30px;
  font-weight: 750;
}
.section {
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 18px;
  padding: 22px;
  margin-top: 18px;
  box-shadow: 0 8px 20px rgba(15, 23, 42, .04);
}
.section h2 {
  margin: 0 0 14px;
  font-size: 21px;
}
.note {
  color: var(--muted);
  font-size: 14px;
}

.legend {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  margin-top: 12px;
}
.legend-item {
  display: flex;
  align-items: center;
  gap: 7px;
  background: #f8fafc;
  border: 1px solid var(--border);
  border-radius: 999px;
  padding: 7px 10px;
  font-size: 13px;
}
.dot {
  width: 11px;
  height: 11px;
  border-radius: 50%;
  display: inline-block;
}
.dot.ok { background: var(--ok); }
.dot.warning { background: var(--warn); }
.dot.problem { background: var(--problem); }
.dot.low { background: var(--low); }
.dot.neutral { background: var(--accent); }

.card-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 14px;
}
.gene-card {
  border: 1px solid var(--border);
  border-radius: 16px;
  background: white;
  padding: 15px;
  box-shadow: 0 6px 16px rgba(15, 23, 42, .04);
}
.gene-card.ok { border-left: 7px solid var(--ok); }
.gene-card.warning { border-left: 7px solid var(--warn); }
.gene-card.problem { border-left: 7px solid var(--problem); }
.gene-card.low { border-left: 7px solid var(--low); }
.gene-card.neutral { border-left: 7px solid var(--accent); }

.gene-card-header {
  display: flex;
  justify-content: space-between;
  gap: 10px;
  align-items: flex-start;
  margin-bottom: 10px;
}
.gene-title {
  font-size: 18px;
  font-weight: 760;
}
.gene-subtitle {
  color: var(--muted);
  font-size: 12px;
  margin-top: 2px;
}
.card-line {
  margin: 7px 0;
  font-size: 13px;
}
.card-line strong {
  color: #334155;
}
.card-text {
  font-size: 13px;
  color: #374151;
  margin-top: 9px;
}
.chips, .source-list {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  margin-top: 8px;
}
.chip {
  background: #f1f5f9;
  color: #334155;
  border-radius: 999px;
  padding: 4px 8px;
  font-size: 12px;
}
.badge {
  display: inline-block;
  border-radius: 999px;
  padding: 3px 9px;
  font-size: 12px;
  font-weight: 650;
  white-space: nowrap;
}
.badge.ok { background: var(--ok-bg); color: var(--ok); }
.badge.warning { background: var(--warn-bg); color: var(--warn); }
.badge.problem { background: var(--problem-bg); color: var(--problem); }
.badge.low { background: var(--low-bg); color: var(--low); }
.badge.neutral { background: var(--accent-bg); color: var(--accent); }
.files {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 12px;
}
.files ul { margin: 8px 0 0; padding-left: 20px; }
code {
  background: #f1f5f9;
  padding: 2px 5px;
  border-radius: 6px;
  font-size: 11px;
}
.footer {
  color: var(--muted);
  margin-top: 22px;
  font-size: 13px;
}
"""

    def badge(value: str) -> str:
        return f'<span class="badge {_status_class(value)}">{_html_escape(value)}</span>'

    def chips(values: str, limit: int = 5) -> str:
        items = [x for x in str(values or ".").split(";") if x and x != "."]
        if not items:
            return ""
        html_items = "".join(f"<span class='chip'>{_html_escape(x)}</span>" for x in items[:limit])
        if len(items) > limit:
            html_items += f"<span class='chip'>+{len(items) - limit} more</span>"
        return f"<div class='chips'>{html_items}</div>"

    def sources(values: str, limit: int = 4) -> str:
        items = [x for x in str(values or ".").split(";") if x and x != "."]
        if not items:
            return ""
        html_items = "".join(f"<code>{_html_escape(x)}</code>" for x in items[:limit])
        if len(items) > limit:
            html_items += f"<code>+{len(items) - limit} more</code>"
        return f"<div class='source-list'>{html_items}</div>"

    def render_inventory_cards(rows: List[Dict[str, str]]) -> str:
        if not rows:
            return "<p class='note'>No gene inventory rows available.</p>"
        cards = ["<div class='card-grid'>"]
        for r in rows:
            status = r.get("annotation_status", ".")
            priority = r.get("annotation_priority", ".")
            card_class = _status_class(status)

            cards.append(
                f"<div class='gene-card {card_class}'>"
                "<div class='gene-card-header'>"
                f"<div><div class='gene-title'>{_html_escape(r.get('gene', '.'))}</div>"
                f"<div class='gene-subtitle'>{_html_escape(r.get('type', '.'))}</div></div>"
                f"<div>{badge(priority)}</div>"
                "</div>"
                f"<div class='card-line'><strong>Status:</strong> {badge(status)}</div>"
                f"<div class='card-line'><strong>Expected:</strong> {badge(r.get('expected_status', '.'))}</div>"
                f"<div class='card-line'><strong>Action:</strong> {_html_escape(r.get('recommended_action', '.'))}</div>"
                f"<div class='card-text'>{_html_escape(r.get('annotation_interpretation', '.'))}</div>"
                f"{chips(r.get('annotation_issues', '.'))}"
                f"{chips(r.get('annotation_recommendations', '.'))}"
                "</div>"
            )
        cards.append("</div>")
        return "\n".join(cards)

    def render_review_cards(rows: List[Dict[str, str]]) -> str:
        if not rows:
            return "<p class='note'>No review targets detected.</p>"
        cards = ["<div class='card-grid'>"]
        for r in rows:
            priority = r.get("annotation_priority", ".")
            card_class = _status_class(priority)

            cards.append(
                f"<div class='gene-card {card_class}'>"
                "<div class='gene-card-header'>"
                f"<div><div class='gene-title'>{_html_escape(r.get('gene', '.'))}</div>"
                f"<div class='gene-subtitle'>{_html_escape(r.get('n_annotation_evidence_sources', '.'))} evidence source(s)</div></div>"
                f"<div>{badge(priority)}</div>"
                "</div>"
                f"<div class='card-text'>{_html_escape(r.get('annotation_interpretation', '.'))}</div>"
                f"{chips(r.get('annotation_issues', '.'))}"
                f"{chips(r.get('annotation_recommendations', '.'))}"
                f"{sources(r.get('annotation_evidence_sources', '.'))}"
                "</div>"
            )
        cards.append("</div>")
        return "\n".join(cards)

    def render_evidence_cards(rows: List[Dict[str, str]], limit: int = 24) -> str:
        if not rows:
            return "<p class='note'>No annotation evidence rows available.</p>"
        cards = ["<div class='card-grid'>"]
        for r in rows[:limit]:
            priority = r.get("priority", ".")
            status = r.get("status", ".")
            card_class = _status_class(priority if priority != "." else status)

            cards.append(
                f"<div class='gene-card {card_class}'>"
                "<div class='gene-card-header'>"
                f"<div><div class='gene-title'>{_html_escape(r.get('gene', '.'))}</div>"
                f"<div class='gene-subtitle'>{_html_escape(r.get('issue_type', '.'))}</div></div>"
                f"<div>{badge(priority)}</div>"
                "</div>"
                f"<div class='card-line'><strong>Status:</strong> {badge(status)}</div>"
                f"<div class='card-line'><strong>Recommendation:</strong> {_html_escape(r.get('recommendation', '.'))}</div>"
                f"<div class='card-line'><strong>Records:</strong> {_html_escape(r.get('n_records', '.'))}</div>"
                f"{sources(r.get('evidence_source', '.'))}"
                "</div>"
            )
        cards.append("</div>")
        if len(rows) > limit:
            cards.append(f"<p class='note'>Showing first {limit} of {len(rows)} evidence items. See TSV for complete output.</p>")
        return "\n".join(cards)

    type_counts = _html_escape(", ".join(f"{k}={v}" for k, v in sorted(by_type.items())) or "not available")
    status_counts = _html_escape(", ".join(f"{k}={v}" for k, v in sorted(by_status.items())) or "not available")

    available_html = "".join(f"<li><code>{_html_escape(x)}</code></li>" for x in available_files) or "<li>None</li>"
    unavailable_html = "".join(f"<li><code>{_html_escape(x)}</code></li>" for x in unavailable_files) or "<li>None</li>"

    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>MitoCurator annotation assessment report</title>
<style>{css}</style>
</head>
<body>
<div class="container">
  <div class="header">
    <h1>MitoCurator annotation assessment report</h1>
    <p>Initial annotation QC, evidence summary and curator-facing recommendations.</p>
  </div>

  <div class="grid">
    <div class="card"><div class="label">Genes/features in inventory</div><div class="value">{total_genes}</div></div>
    <div class="card"><div class="label">OK / keep</div><div class="value">{ok_genes}</div></div>
    <div class="card"><div class="label">Problematic</div><div class="value">{problem_genes}</div></div>
    <div class="card"><div class="label">Missing/incomplete</div><div class="value">{missing_genes}</div></div>
    <div class="card"><div class="label">Review targets</div><div class="value">{review_targets}</div></div>
    <div class="card"><div class="label">High-priority targets</div><div class="value">{high_targets}</div></div>
  </div>

  <div class="section">
    <h2>Legend and overview</h2>
    <p class="note">Type counts: {type_counts}</p>
    <p class="note">Status counts: {status_counts}</p>
    <div class="legend">
      <div class="legend-item"><span class="dot ok"></span> OK / keep</div>
      <div class="legend-item"><span class="dot problem"></span> High-priority, missing or problematic</div>
      <div class="legend-item"><span class="dot warning"></span> Review / candidate / medium priority</div>
      <div class="legend-item"><span class="dot low"></span> Low-priority signal</div>
      <div class="legend-item"><span class="dot neutral"></span> Informational / not classified</div>
    </div>
  </div>

  <div class="section">
    <h2>Review targets</h2>
    <p class="note">Genes/features requiring curator attention before downstream read-backed correction.</p>
    {render_review_cards(review_rows)}
  </div>

  <div class="section">
    <h2>Complete annotation gene inventory</h2>
    <p class="note">All expected genes/features, including OK/KEEP entries.</p>
    {render_inventory_cards(inventory_rows)}
  </div>

  <div class="section">
    <h2>Annotation evidence summary</h2>
    <p class="note">Aggregated annotation-level evidence used to classify genes/features. The complete machine-readable table is available as TSV.</p>
    {render_evidence_cards(evidence_rows)}
  </div>

  <div class="section">
    <h2>Evidence files</h2>
    <div class="files">
      <div>
        <h3>Available</h3>
        <ul>{available_html}</ul>
      </div>
      <div>
        <h3>Not available</h3>
        <ul>{unavailable_html}</ul>
      </div>
    </div>
  </div>

  <div class="section">
    <h2>Output files</h2>
    <ul>
      <li><code>{_html_escape(md_path.name)}</code> — Markdown report</li>
      <li><code>{_html_escape(annotation_summary_tsv.name)}</code> — annotation evidence summary</li>
      <li><code>{_html_escape(gene_inventory_tsv.name)}</code> — complete gene inventory</li>
      <li><code>{_html_escape(gene_summary_tsv.name)}</code> — review targets</li>
    </ul>
  </div>

  <div class="footer">
    Diagnostic-only report. No GenBank file is modified by this assessment.
  </div>
</div>
</body>
</html>
"""
    html_path.write_text(html_text, encoding="utf-8")


def _write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def generate_annotation_assessment_report(root: Path):
    """Generate an annotation/refinement report.

    This report intentionally uses only annotation-level evidence produced before
    downstream read-backed curation stages. It does not summarize read support,
    targeted consensus, reconstruction pools or candidate assembly.
    """
    report_dir = ensure_dir(root / "06_annotation_assessment")
    md_path = report_dir / "annotation_assessment_report.md"
    html_path = report_dir / "annotation_assessment_report.html"
    annotation_summary_tsv = report_dir / "annotation_evidence_summary.tsv"
    gene_summary_tsv = report_dir / "annotation_review_targets.tsv"
    gene_inventory_tsv = report_dir / "annotation_gene_inventory.tsv"

    evidence = _collect_evidence_status(root)
    raw_rows = _collect_annotation_rows(root)
    collapsed_rows = _collapse_annotation_rows(raw_rows)
    inventory_rows = _build_gene_inventory(collapsed_rows)
    gene_rows = _build_gene_annotation_summary(collapsed_rows)

    _write_tsv(
        annotation_summary_tsv,
        collapsed_rows,
        [
            "gene",
            "issue_type",
            "status",
            "priority",
            "recommendation",
            "evidence_source",
            "n_records",
            "evidence_examples",
        ],
    )

    _write_tsv(
        gene_inventory_tsv,
        inventory_rows,
        [
            "gene",
            "type",
            "expected_status",
            "annotation_status",
            "annotation_priority",
            "recommended_action",
            "annotation_issues",
            "annotation_recommendations",
            "annotation_evidence_sources",
            "annotation_interpretation",
        ],
    )

    _write_tsv(
        gene_summary_tsv,
        gene_rows,
        [
            "gene",
            "annotation_priority",
            "n_annotation_evidence_sources",
            "annotation_evidence_sources",
            "annotation_issues",
            "annotation_recommendations",
            "annotation_interpretation",
        ],
    )

    available = [r["file"] for r in evidence if r["status"] == "available"]
    unavailable = [r["file"] for r in evidence if r["status"] != "available"]

    with open(md_path, "w", encoding="utf-8") as out:
        out.write("# MitoCurator annotation/refinement report\n\n")

        out.write("## 1. Scope\n\n")
        out.write(
            "This report summarizes annotation-level evidence produced before downstream "
            "read-backed curation stages. It is intended to help the researcher decide which "
            "genes or features require manual inspection before read-support, targeted "
            "consensus or local assembly steps are interpreted.\n\n"
        )
        out.write("This report does **not** summarize read support, targeted consensus, reconstruction pools or candidate assembly.\n\n")

        out.write("## 2. Annotation evidence files\n\n")
        out.write("### Available\n\n")
        if available:
            for item in available:
                out.write(f"- `{item}`\n")
        else:
            out.write("- No annotation evidence files available.\n")

        out.write("\n### Not available\n\n")
        if unavailable:
            for item in unavailable:
                out.write(f"- `{item}`\n")
        else:
            out.write("- All expected annotation evidence files are available.\n")

        out.write("\n## 3. Annotation/refinement summary\n\n")
        out.write(f"- Raw annotation evidence rows: `{len(raw_rows)}`\n")
        out.write(f"- Aggregated annotation evidence rows: `{len(collapsed_rows)}`\n")
        out.write(f"- Genes/features in complete inventory: `{len(inventory_rows)}`\n")
        out.write(f"- Genes/features requiring review: `{len(gene_rows)}`\n")
        out.write(f"- Evidence summary table: `{annotation_summary_tsv.name}`\n")
        out.write(f"- Complete gene inventory: `{gene_inventory_tsv.name}`\n")
        out.write(f"- Review targets table: `{gene_summary_tsv.name}`\n\n")

        if not collapsed_rows:
            out.write("No annotation/refinement issues were summarized from the available inputs.\n")
        else:
            out.write("| gene | issue | priority | recommendation | evidence source | n |\n")
            out.write("|---|---|---|---|---|---:|\n")
            for r in collapsed_rows:
                out.write(
                    f"| {r['gene']} | {r['issue_type']} | {r['priority']} | "
                    f"{r['recommendation']} | `{r['evidence_source']}` | {r['n_records']} |\n"
                )

        out.write("\n## 4. Complete annotation gene inventory\n\n")
        out.write(f"Complete table: `{gene_inventory_tsv.name}`\n\n")
        if not inventory_rows:
            out.write("No gene inventory was generated.\n")
        else:
            out.write("| gene | type | status | priority | action |\n")
            out.write("|---|---|---|---|---|\n")
            for r in inventory_rows:
                out.write(
                    f"| {r['gene']} | {r['type']} | {r['annotation_status']} | "
                    f"{r['annotation_priority']} | {r['recommended_action']} |\n"
                )

        out.write("\n## 5. Review targets and recommended next steps\n\n")
        if not gene_rows:
            out.write("No gene-level annotation review targets were generated.\n")
        else:
            out.write("| gene | priority | evidence sources | interpretation |\n")
            out.write("|---|---|---:|---|\n")
            for r in gene_rows:
                out.write(
                    f"| {r['gene']} | {r['annotation_priority']} | "
                    f"{r['n_annotation_evidence_sources']} | {r['annotation_interpretation']} |\n"
                )

        out.write("\n## 6. Diagnostic-only statement\n\n")
        out.write(
            "This report is **diagnostic-only**. It does not alter the GenBank file, "
            "does not accept or reject candidate corrections automatically, and does not replace manual curation.\n"
        )

    _write_html_report(
        html_path=html_path,
        evidence=evidence,
        inventory_rows=inventory_rows,
        review_rows=gene_rows,
        evidence_rows=collapsed_rows,
        md_path=md_path,
        annotation_summary_tsv=annotation_summary_tsv,
        gene_inventory_tsv=gene_inventory_tsv,
        gene_summary_tsv=gene_summary_tsv,
    )

    return md_path, annotation_summary_tsv


# Backwards-compatible alias for older internal calls.
def generate_final_report(root: Path):
    return generate_annotation_assessment_report(root)
