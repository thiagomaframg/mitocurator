from __future__ import annotations

import csv
import html
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .utils import ensure_dir


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, "r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _escape(value: object) -> str:
    return html.escape(str(value if value is not None else "."))


def _first(rows: List[Dict[str, str]], key: str, default: str = ".") -> str:
    for row in rows:
        value = str(row.get(key, "") or "")
        if value and value != ".":
            return value
    return default


def _join_unique(values: List[str], sep: str = ";") -> str:
    seen = []
    for value in values:
        value = str(value or ".")
        if not value or value == ".":
            continue
        for part in value.split(";"):
            part = part.strip()
            if part and part != "." and part not in seen:
                seen.append(part)
    return sep.join(seen) if seen else "."


def _num(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _index_rows(paths: Dict[str, Path]) -> List[Dict[str, str]]:
    rows = []
    for label, path in paths.items():
        rows.append({
            "evidence_label": label,
            "path": str(path),
            "status": "available" if path.exists() and path.stat().st_size > 0 else "not_available",
            "size_bytes": str(path.stat().st_size) if path.exists() else "0",
        })
    return rows


def _group_by_gene(rows: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
    grouped: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        gene = str(row.get("gene", ".") or ".")
        grouped.setdefault(gene, []).append(row)
    return grouped


def _variant_summary_for_gene(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return "."
    parts = []
    for row in rows:
        parts.append(
            f"{row.get('readset', '.')}:"
            f"{row.get('n_variants', '0')} variants/"
            f"{row.get('n_snps', '0')} SNPs/"
            f"{row.get('n_indels', '0')} indels/"
            f"maxAF={row.get('max_alt_frequency', '.')}"
        )
    return "; ".join(parts)


def _coverage_summary_for_gene(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return "."
    parts = []
    for row in rows:
        parts.append(
            f"{row.get('readset', '.')}:"
            f"mean={row.get('mean_depth', '.')},"
            f"min={row.get('min_depth', '.')},"
            f"covered={row.get('pct_bases_covered', '.')}%"
        )
    return "; ".join(parts)


def _read_support_summary_for_gene(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return "."
    parts = []
    for row in rows:
        parts.append(
            f"stop{row.get('stop_aa_position', '.')}:"
            f"{row.get('consensus_recommendation', '.')}"
            f"({row.get('priority', '.')})"
        )
    return "; ".join(parts)


def _targeted_consensus_summary_for_gene(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return "."
    parts = []
    for row in rows[:5]:
        parts.append(
            f"{row.get('target_id', '.')}:"
            f"{row.get('combined_recommendation', row.get('recommendation', '.'))},"
            f"rank={row.get('combined_rank', row.get('rank', '.'))},"
            f"readsets={row.get('read_sets_supporting', row.get('read_set', '.'))}"
        )
    return "; ".join(parts)


def _decide_gene(
    gene: str,
    inventory: List[Dict[str, str]],
    review: List[Dict[str, str]],
    coverage: List[Dict[str, str]],
    variants: List[Dict[str, str]],
    read_support: List[Dict[str, str]],
    cross_consensus: List[Dict[str, str]],
    best_consensus: List[Dict[str, str]],
) -> Tuple[str, str, str, str]:
    issues = _join_unique([r.get("annotation_issues", ".") for r in inventory + review])
    annotation_status = _first(inventory, "annotation_status", ".")
    annotation_priority = _first(review, "annotation_priority", _first(inventory, "annotation_priority", "."))

    rs_recs = {r.get("consensus_recommendation", ".") for r in read_support}
    cross_recs = {r.get("combined_recommendation", ".") for r in cross_consensus}
    best_recs = {r.get("recommendation", ".") for r in best_consensus}

    has_good_cross_gene_candidate = any("GENE_CANDIDATE_SUPPORTED" in x for x in cross_recs)
    has_good_best_gene_candidate = any("GENE_CANDIDATE_SUPPORTED" in x or "PARTIAL_GENE_CANDIDATE_SUPPORTED" in x for x in best_recs)

    n_indels = sum(int(_num(r.get("n_indels", "0"))) for r in variants)
    max_af = max([_num(r.get("max_alt_frequency", "0")) for r in variants] or [0.0])
    min_cov = min([_num(r.get("min_depth", "0")) for r in coverage] or [0.0])
    mean_cov = max([_num(r.get("mean_depth", "0")) for r in coverage] or [0.0])

    if has_good_cross_gene_candidate:
        return (
            "MISSING_GENE_CANDIDATE_SUPPORTED_ACROSS_READSETS",
            "HIGH",
            "Add or review candidate gene annotation using cross-readset targeted consensus evidence.",
            "Missing/partial gene candidate is supported by targeted consensus in multiple readsets.",
        )

    if has_good_best_gene_candidate:
        return (
            "MISSING_GENE_CANDIDATE_SUPPORTED_BY_TARGETED_CONSENSUS",
            "HIGH",
            "Review candidate gene coordinates and consensus FASTA before adding annotation.",
            "Missing/partial gene candidate is supported by targeted consensus in at least one readset.",
        )

    if "CORRECTION_SUPPORTED_BY_ALL_READSETS" in rs_recs:
        return (
            "SEQUENCE_CORRECTION_SUPPORTED_BY_ALL_READSETS",
            "HIGH",
            "Review local sequence correction before updating GenBank/FASTA.",
            "All evaluated readsets support replacing a problematic stop codon with a coding codon.",
        )

    if "CORRECTION_SUPPORTED_BY_SOME_READSETS" in rs_recs or "CONFLICTING_READSET_EVIDENCE" in rs_recs:
        return (
            "SEQUENCE_CORRECTION_REQUIRES_MANUAL_REVIEW",
            "MEDIUM",
            "Inspect read-level evidence, coverage and strand/readset consistency.",
            "Readsets provide partial or conflicting support for a local correction.",
        )

    if "STOP_CONFIRMED_BY_ALL_READSETS" in rs_recs:
        return (
            "STOP_CONFIRMED_BY_READS",
            "HIGH",
            "Treat as likely biological/annotation-boundary issue; inspect CDS boundaries and reference alignment.",
            "All evaluated readsets support the stop codon present in the current sequence.",
        )

    if "STOP_SUPPORTED_BY_SOME_READSETS" in rs_recs:
        return (
            "STOP_PARTIALLY_SUPPORTED_BY_READS",
            "MEDIUM",
            "Inspect readset-specific evidence and consider manual review.",
            "At least one readset supports the stop codon, but evidence is not unanimous.",
        )

    if annotation_status == "MISSING_OR_INCOMPLETE":
        return (
            "MISSING_OR_INCOMPLETE_UNRESOLVED",
            "HIGH",
            "Prioritize targeted consensus/local assembly or manual annotation review.",
            "Annotation assessment reports missing or incomplete feature without decisive downstream support.",
        )

    if annotation_status == "PROBLEMATIC":
        if n_indels > 0 and max_af >= 0.8:
            return (
                "PROBLEMATIC_WITH_HIGH_FREQUENCY_INDEL",
                "HIGH",
                "Inspect variant evidence for possible frameshift/assembly correction.",
                "Problematic feature overlaps high-frequency indel evidence.",
            )
        return (
            "PROBLEMATIC_REQUIRES_MANUAL_REVIEW",
            annotation_priority if annotation_priority != "." else "MEDIUM",
            "Review annotation boundaries, translation frame, reference similarity and read evidence.",
            "Annotation-level problem remains after integrated evidence review.",
        )

    if annotation_status == "OK" and mean_cov > 0 and min_cov > 0:
        return (
            "NO_ACTION_SUPPORTED",
            "LOW",
            "No immediate curation action required.",
            "Feature is annotated as OK and has read coverage support.",
        )

    return (
        "INSUFFICIENT_EVIDENCE",
        "MEDIUM",
        "Review available evidence manually.",
        "Integrated report could not assign a decisive interpretation.",
    )


def _build_gene_decisions(root: Path) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    paths = {
        "annotation_review_targets": root / "06_annotation_assessment" / "annotation_review_targets.tsv",
        "annotation_gene_inventory": root / "06_annotation_assessment" / "annotation_gene_inventory.tsv",
        "mapping_summary": root / "08_read_mapping" / "mapping_summary.tsv",
        "coverage_by_gene": root / "08_read_mapping" / "coverage_by_gene.tsv",
        "gene_variant_evidence": root / "09_variant_evidence" / "gene_variant_evidence.tsv",
        "readset_consensus_recommendations": root / "10_read_support" / "readset_consensus_recommendations.tsv",
        "best_missing_gene_candidates": root / "13_targeted_consensus" / "best_missing_gene_candidates.tsv",
        "cross_readset_missing_gene_candidates": root / "13_targeted_consensus" / "cross_readset_missing_gene_candidates.tsv",
    }

    evidence_index = _index_rows(paths)

    review = _read_tsv(paths["annotation_review_targets"])
    inventory = _read_tsv(paths["annotation_gene_inventory"])
    coverage = _read_tsv(paths["coverage_by_gene"])
    variants = _read_tsv(paths["gene_variant_evidence"])
    read_support = _read_tsv(paths["readset_consensus_recommendations"])
    best_consensus = _read_tsv(paths["best_missing_gene_candidates"])
    cross_consensus = _read_tsv(paths["cross_readset_missing_gene_candidates"])

    by_inv = _group_by_gene(inventory)
    by_review = _group_by_gene(review)
    by_cov = _group_by_gene(coverage)
    by_var = _group_by_gene(variants)
    by_rs = _group_by_gene(read_support)
    by_best = _group_by_gene(best_consensus)
    by_cross = _group_by_gene(cross_consensus)

    genes = sorted(
        set(by_inv)
        | set(by_review)
        | set(by_cov)
        | set(by_var)
        | set(by_rs)
        | set(by_best)
        | set(by_cross)
    )

    rows: List[Dict[str, str]] = []

    for gene in genes:
        inv = by_inv.get(gene, [])
        rev = by_review.get(gene, [])
        cov = by_cov.get(gene, [])
        var = by_var.get(gene, [])
        rs = by_rs.get(gene, [])
        best = by_best.get(gene, [])
        cross = by_cross.get(gene, [])

        decision, priority, action, interpretation = _decide_gene(
            gene, inv, rev, cov, var, rs, cross, best
        )

        rows.append({
            "gene": gene,
            "feature_type": _first(inv, "type", _first(cov, "type", _first(var, "feature_type", "."))),
            "initial_annotation_status": _first(inv, "annotation_status", "."),
            "initial_annotation_priority": _first(rev, "annotation_priority", _first(inv, "annotation_priority", ".")),
            "integrated_decision": decision,
            "integrated_priority": priority,
            "recommended_action": action,
            "integrated_interpretation": interpretation,
            "annotation_issues": _join_unique([r.get("annotation_issues", ".") for r in inv + rev]),
            "annotation_recommendations": _join_unique([r.get("annotation_recommendations", ".") for r in inv + rev]),
            "coverage_summary": _coverage_summary_for_gene(cov),
            "variant_summary": _variant_summary_for_gene(var),
            "read_support_summary": _read_support_summary_for_gene(rs),
            "targeted_consensus_summary": _targeted_consensus_summary_for_gene(cross or best),
        })

    return rows, evidence_index


def _status_class(value: str) -> str:
    v = str(value or "").upper()
    if "SUPPORTED" in v or "NO_ACTION" in v:
        return "ok"
    if "HIGH" in v or "CORRECTION" in v or "MISSING" in v or "PROBLEMATIC" in v:
        return "problem"
    if "MEDIUM" in v or "REVIEW" in v or "PARTIAL" in v or "INSUFFICIENT" in v:
        return "warning"
    if "LOW" in v:
        return "low"
    return "neutral"


def _write_markdown(path: Path, rows: List[Dict[str, str]], evidence_index: List[Dict[str, str]]) -> None:
    n_high = sum(1 for r in rows if r["integrated_priority"] == "HIGH")
    n_med = sum(1 for r in rows if r["integrated_priority"] == "MEDIUM")
    n_low = sum(1 for r in rows if r["integrated_priority"] == "LOW")

    lines = []
    lines.append("# MitoCurator integrated curation report\n")
    lines.append("This report integrates annotation assessment, read mapping, variant evidence, read-support interpretation and targeted consensus evidence.\n")
    lines.append("## Summary\n")
    lines.append(f"- Genes/features summarized: {len(rows)}")
    lines.append(f"- High-priority decisions: {n_high}")
    lines.append(f"- Medium-priority decisions: {n_med}")
    lines.append(f"- Low-priority decisions: {n_low}\n")

    lines.append("## Evidence index\n")
    lines.append("| Evidence | Status | Path |")
    lines.append("|---|---:|---|")
    for e in evidence_index:
        lines.append(f"| {e['evidence_label']} | {e['status']} | `{e['path']}` |")
    lines.append("")

    lines.append("## Integrated gene decisions\n")
    lines.append("| Gene | Decision | Priority | Recommended action |")
    lines.append("|---|---|---:|---|")
    for r in rows:
        lines.append(
            f"| {r['gene']} | {r['integrated_decision']} | {r['integrated_priority']} | {r['recommended_action']} |"
        )
    lines.append("")

    lines.append("## Detailed evidence by gene\n")
    for r in rows:
        lines.append(f"### {r['gene']}\n")
        lines.append(f"- Feature type: {r['feature_type']}")
        lines.append(f"- Initial annotation status: {r['initial_annotation_status']}")
        lines.append(f"- Integrated decision: **{r['integrated_decision']}**")
        lines.append(f"- Priority: **{r['integrated_priority']}**")
        lines.append(f"- Recommended action: {r['recommended_action']}")
        lines.append(f"- Interpretation: {r['integrated_interpretation']}")
        lines.append(f"- Annotation issues: {r['annotation_issues']}")
        lines.append(f"- Coverage: {r['coverage_summary']}")
        lines.append(f"- Variant evidence: {r['variant_summary']}")
        lines.append(f"- Read support: {r['read_support_summary']}")
        lines.append(f"- Targeted consensus: {r['targeted_consensus_summary']}\n")

    path.write_text("\n".join(lines), encoding="utf-8")


def _write_html(path: Path, rows: List[Dict[str, str]], evidence_index: List[Dict[str, str]]) -> None:
    n_high = sum(1 for r in rows if r["integrated_priority"] == "HIGH")
    n_med = sum(1 for r in rows if r["integrated_priority"] == "MEDIUM")
    n_low = sum(1 for r in rows if r["integrated_priority"] == "LOW")

    css = """
body { margin:0; background:#f6f8fb; color:#1f2937; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Arial,sans-serif; }
.container { max-width:1200px; margin:0 auto; padding:32px 20px 56px; }
.header { background:linear-gradient(135deg,#065f46,#2563eb); color:white; border-radius:22px; padding:28px 32px; box-shadow:0 12px 30px rgba(15,23,42,.16); }
.header h1 { margin:0 0 8px; font-size:30px; }
.header p { margin:0; opacity:.92; }
.grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); gap:14px; margin:20px 0; }
.card, .section, .gene-card { background:white; border:1px solid #e5e7eb; border-radius:18px; box-shadow:0 8px 20px rgba(15,23,42,.05); }
.card { padding:18px; }
.card .label { color:#6b7280; font-size:13px; margin-bottom:8px; }
.card .value { font-size:30px; font-weight:750; }
.section { padding:22px; margin-top:18px; }
.section h2 { margin:0 0 14px; font-size:21px; }
.note { color:#6b7280; font-size:14px; }
.legend { display:flex; flex-wrap:wrap; gap:10px; margin-top:12px; }
.legend-item { display:flex; align-items:center; gap:7px; background:#f8fafc; border:1px solid #e5e7eb; border-radius:999px; padding:7px 10px; font-size:13px; }
.dot { width:11px; height:11px; border-radius:50%; display:inline-block; }
.dot.ok { background:#047857; } .dot.warning { background:#b45309; } .dot.problem { background:#b91c1c; } .dot.low { background:#475569; } .dot.neutral { background:#1d4ed8; }
.gene-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); gap:14px; }
.gene-card { padding:15px; border-left:7px solid #1d4ed8; }
.gene-card.ok { border-left-color:#047857; }
.gene-card.warning { border-left-color:#b45309; }
.gene-card.problem { border-left-color:#b91c1c; }
.gene-card.low { border-left-color:#475569; }
.gene-title { font-size:18px; font-weight:760; margin-bottom:4px; }
.badge { display:inline-block; border-radius:999px; padding:3px 9px; font-size:12px; font-weight:650; }
.badge.ok { background:#d1fae5; color:#047857; }
.badge.warning { background:#fef3c7; color:#b45309; }
.badge.problem { background:#fee2e2; color:#b91c1c; }
.badge.low { background:#e2e8f0; color:#475569; }
.badge.neutral { background:#dbeafe; color:#1d4ed8; }
.line { margin:7px 0; font-size:13px; }
.line strong { color:#334155; }
.small { color:#475569; font-size:12px; overflow-wrap:anywhere; }
code { background:#f1f5f9; padding:2px 5px; border-radius:6px; font-size:11px; }
table { width:100%; border-collapse:collapse; font-size:13px; }
th, td { padding:8px; border-bottom:1px solid #e5e7eb; text-align:left; vertical-align:top; }
"""

    def badge(value: str) -> str:
        cls = _status_class(value)
        return f'<span class="badge {cls}">{_escape(value)}</span>'

    cards = []
    for r in rows:
        cls = _status_class(r["integrated_decision"] + " " + r["integrated_priority"])
        cards.append(f"""
<div class="gene-card {cls}">
  <div class="gene-title">{_escape(r['gene'])} {badge(r['integrated_priority'])}</div>
  <div class="line"><strong>Decision:</strong> {badge(r['integrated_decision'])}</div>
  <div class="line"><strong>Action:</strong> {_escape(r['recommended_action'])}</div>
  <div class="line"><strong>Initial status:</strong> {_escape(r['initial_annotation_status'])}</div>
  <div class="line"><strong>Interpretation:</strong> {_escape(r['integrated_interpretation'])}</div>
  <div class="small"><strong>Coverage:</strong> {_escape(r['coverage_summary'])}</div>
  <div class="small"><strong>Variants:</strong> {_escape(r['variant_summary'])}</div>
  <div class="small"><strong>Read support:</strong> {_escape(r['read_support_summary'])}</div>
  <div class="small"><strong>Targeted consensus:</strong> {_escape(r['targeted_consensus_summary'])}</div>
</div>
""")

    evidence_rows = "\n".join(
        f"<tr><td>{_escape(e['evidence_label'])}</td><td>{badge(e['status'])}</td><td><code>{_escape(e['path'])}</code></td></tr>"
        for e in evidence_index
    )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>MitoCurator integrated curation report</title>
<style>{css}</style>
</head>
<body>
<div class="container">
  <div class="header">
    <h1>MitoCurator integrated curation report</h1>
    <p>Annotation, mapping, variant, read-support and targeted-consensus evidence integrated into gene-level decisions.</p>
  </div>

  <div class="grid">
    <div class="card"><div class="label">Genes/features summarized</div><div class="value">{len(rows)}</div></div>
    <div class="card"><div class="label">High priority</div><div class="value">{n_high}</div></div>
    <div class="card"><div class="label">Medium priority</div><div class="value">{n_med}</div></div>
    <div class="card"><div class="label">Low priority</div><div class="value">{n_low}</div></div>
  </div>

  <div class="section">
    <h2>Legend</h2>
    <div class="legend">
      <div class="legend-item"><span class="dot ok"></span> Supported/no-action evidence</div>
      <div class="legend-item"><span class="dot problem"></span> High-priority correction or missing/probl. feature</div>
      <div class="legend-item"><span class="dot warning"></span> Manual review / unresolved</div>
      <div class="legend-item"><span class="dot low"></span> Low-priority evidence</div>
      <div class="legend-item"><span class="dot neutral"></span> Informational</div>
    </div>
  </div>

  <div class="section">
    <h2>Integrated gene decisions</h2>
    <div class="gene-grid">
      {''.join(cards)}
    </div>
  </div>

  <div class="section">
    <h2>Evidence index</h2>
    <table>
      <thead><tr><th>Evidence</th><th>Status</th><th>Path</th></tr></thead>
      <tbody>{evidence_rows}</tbody>
    </table>
  </div>
</div>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def run_integrated_report(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "15_integrated_report")

    decision_rows, evidence_index = _build_gene_decisions(root)

    decision_fields = [
        "gene",
        "feature_type",
        "initial_annotation_status",
        "initial_annotation_priority",
        "integrated_decision",
        "integrated_priority",
        "recommended_action",
        "integrated_interpretation",
        "annotation_issues",
        "annotation_recommendations",
        "coverage_summary",
        "variant_summary",
        "read_support_summary",
        "targeted_consensus_summary",
    ]

    _write_tsv(outdir / "integrated_gene_decisions.tsv", decision_rows, decision_fields)
    _write_tsv(outdir / "evidence_index.tsv", evidence_index, ["evidence_label", "path", "status", "size_bytes"])
    _write_markdown(outdir / "integrated_curation_report.md", decision_rows, evidence_index)
    _write_html(outdir / "integrated_curation_report.html", decision_rows, evidence_index)

    return outdir
