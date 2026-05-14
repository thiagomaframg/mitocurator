from __future__ import annotations

import csv
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
    annotation_summary_tsv = report_dir / "annotation_evidence_summary.tsv"
    gene_summary_tsv = report_dir / "annotation_review_targets.tsv"

    evidence = _collect_evidence_status(root)
    raw_rows = _collect_annotation_rows(root)
    collapsed_rows = _collapse_annotation_rows(raw_rows)
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
        out.write(f"- Genes/features with annotation signals: `{len(gene_rows)}`\n")
        out.write(f"- Detailed table: `{annotation_summary_tsv.name}`\n")
        out.write(f"- Gene-level table: `{gene_summary_tsv.name}`\n\n")

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

        out.write("\n## 4. Gene-level annotation interpretation\n\n")
        if not gene_rows:
            out.write("No gene-level annotation interpretation was generated.\n")
        else:
            out.write("| gene | priority | evidence sources | interpretation |\n")
            out.write("|---|---|---:|---|\n")
            for r in gene_rows:
                out.write(
                    f"| {r['gene']} | {r['annotation_priority']} | "
                    f"{r['n_annotation_evidence_sources']} | {r['annotation_interpretation']} |\n"
                )

        out.write("\n## 5. Diagnostic-only statement\n\n")
        out.write(
            "This report is **diagnostic-only**. It does not alter the GenBank file, "
            "does not accept or reject candidate corrections automatically, and does not replace manual curation.\n"
        )

    return md_path, annotation_summary_tsv


# Backwards-compatible alias for older internal calls.
def generate_final_report(root: Path):
    return generate_annotation_assessment_report(root)
