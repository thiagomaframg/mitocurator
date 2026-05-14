from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

from .utils import ensure_dir


EVIDENCE_PATHS = [
    "05_refinement/curation_recommendations.tsv",
    "05_refinement/curation_recommendations.md",
    "06_read_support/read_support_summary.md",
    "06_read_support/readset_consensus_recommendations.tsv",
    "06_read_support/readset_consensus_recommendations.md",
    "10_targeted_consensus/targeted_consensus.tsv",
    "10_targeted_consensus/best_missing_gene_candidates.tsv",
    "10_targeted_consensus/cross_readset_missing_gene_candidates.tsv",
    "11_candidate_assembly/candidate_assembly_summary.tsv",
    "11_candidate_assembly/candidate_assembly_summary.md",
    "07_gene_qc/gene_qc.tsv",
    "07_gene_qc/problematic_features.tsv",
    "07_gene_qc/diagnostic_summary.md",
]


def _first_existing_column(fieldnames: List[str], names: List[str]) -> str | None:
    for name in names:
        if name in fieldnames:
            return name
    return None


def _safe_read_tsv(path: Path) -> List[Dict[str, str]]:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            return list(reader)
    except Exception:
        return []


def _collect_evidence_status(root: Path) -> List[Dict[str, str]]:
    rows = []
    for rel in EVIDENCE_PATHS:
        p = root / rel
        rows.append({"file": rel, "status": "available" if p.exists() else "not available"})
    cand_glob = list((root / "11_candidate_assembly").glob("*/diagnosis/candidate_gene_diagnosis.tsv"))
    if cand_glob:
        for p in sorted(cand_glob):
            rows.append({"file": str(p.relative_to(root)), "status": "available"})
    else:
        rows.append({"file": "11_candidate_assembly/*/diagnosis/candidate_gene_diagnosis.tsv", "status": "not available"})
    return rows


def _build_gene_problem_priority(root: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []

    p_problematic = root / "07_gene_qc/problematic_features.tsv"
    if p_problematic.exists():
        data = _safe_read_tsv(p_problematic)
        if data:
            cols = list(data[0].keys())
            gene_col = _first_existing_column(cols, ["gene_normalized", "gene", "feature", "feature_name"])
            hint_col = _first_existing_column(cols, ["decision_hint", "status", "comment"])
            if gene_col:
                for r in data:
                    rows.append({"gene": r.get(gene_col, "."), "problem": "problematic_feature", "priority": "high", "evidence_source": "07_gene_qc/problematic_features.tsv", "detail": r.get(hint_col, ".") if hint_col else "."})

    p_refine = root / "05_refinement/curation_recommendations.tsv"
    if p_refine.exists():
        data = _safe_read_tsv(p_refine)
        if data:
            cols = list(data[0].keys())
            gene_col = _first_existing_column(cols, ["gene", "gene_normalized", "target_gene"])
            problem_col = _first_existing_column(cols, ["problem", "issue", "category", "recommendation_type"])
            pri_col = _first_existing_column(cols, ["priority", "decision_hint", "status"])
            detail_col = _first_existing_column(cols, ["recommendation", "comment", "details"])
            for r in data:
                rows.append({"gene": r.get(gene_col, ".") if gene_col else ".", "problem": r.get(problem_col, "refinement_recommendation") if problem_col else "refinement_recommendation", "priority": r.get(pri_col, "review") if pri_col else "review", "evidence_source": "05_refinement/curation_recommendations.tsv", "detail": r.get(detail_col, ".") if detail_col else "."})

    p_targeted = root / "10_targeted_consensus/targeted_consensus.tsv"
    if p_targeted.exists():
        data = _safe_read_tsv(p_targeted)
        if data:
            cols = list(data[0].keys())
            gene_col = _first_existing_column(cols, ["gene", "target_gene", "gene_normalized"])
            status_col = _first_existing_column(cols, ["status", "consensus_status", "decision_hint"])
            for r in data:
                rows.append({"gene": r.get(gene_col, ".") if gene_col else ".", "problem": "targeted_consensus_signal", "priority": r.get(status_col, "review") if status_col else "review", "evidence_source": "10_targeted_consensus/targeted_consensus.tsv", "detail": r.get(status_col, ".") if status_col else "."})

    p_cand = root / "11_candidate_assembly/candidate_assembly_summary.tsv"
    if p_cand.exists():
        data = _safe_read_tsv(p_cand)
        if data:
            cols = list(data[0].keys())
            gene_col = _first_existing_column(cols, ["gene", "target_gene", "gene_normalized"])
            status_col = _first_existing_column(cols, ["status", "decision_hint", "priority"])
            detail_col = _first_existing_column(cols, ["comment", "recommendation", "summary"])
            for r in data:
                rows.append({"gene": r.get(gene_col, ".") if gene_col else ".", "problem": "candidate_assembly_signal", "priority": r.get(status_col, "review") if status_col else "review", "evidence_source": "11_candidate_assembly/candidate_assembly_summary.tsv", "detail": r.get(detail_col, ".") if detail_col else "."})

    return rows


def generate_final_report(root: Path):
    final_dir = ensure_dir(root / "12_final_report")
    md_path = final_dir / "final_curation_report.md"
    tsv_path = final_dir / "final_curation_summary.tsv"

    evidence = _collect_evidence_status(root)
    summary_rows = _build_gene_problem_priority(root)

    with open(tsv_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["gene", "problem", "priority", "evidence_source", "detail"], delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    available = [r["file"] for r in evidence if r["status"] == "available"]
    unavailable = [r["file"] for r in evidence if r["status"] != "available"]

    with open(md_path, "w", encoding="utf-8") as out:
        out.write("# MitoCurator final curation report\n\n")
        out.write("## 1. Visão geral da análise\n\n")
        out.write("Este relatório agrega evidências já produzidas pelas etapas anteriores sem alterar algoritmos ou reanotar automaticamente o GenBank.\n\n")
        out.write("## 2. Arquivos de evidência encontrados\n\n")
        for item in available:
            out.write(f"- `{item}`: available\n")
        for item in unavailable:
            out.write(f"- `{item}`: not available\n")

        out.write("\n## 3. Resumo por gene/problema/prioridade\n\n")
        if not summary_rows:
            out.write("Nenhum resumo tabular disponível; arquivos de entrada podem estar ausentes neste run.\n")
        else:
            grouped: Dict[tuple, int] = {}
            for row in summary_rows:
                key = (row["gene"], row["problem"], row["priority"])
                grouped[key] = grouped.get(key, 0) + 1
            out.write("| gene | problem | priority | count |\n|---|---|---|---:|\n")
            for (gene, problem, priority), count in sorted(grouped.items()):
                out.write(f"| {gene} | {problem} | {priority} | {count} |\n")

        out.write("\n## 4. Recomendações finais integradas\n\n")
        if not summary_rows:
            out.write("- Not available: sem sinais consolidados de refinement/read support/targeted consensus/candidate assembly neste diretório.\n")
        else:
            src_counts: Dict[str, int] = {}
            for row in summary_rows:
                src = row["evidence_source"]
                src_counts[src] = src_counts.get(src, 0) + 1
            out.write("- Priorizar revisão manual de genes com múltiplos sinais em fontes independentes.\n")
            out.write("- Consolidar decisões considerando, quando disponível: refinement, read support, targeted consensus e candidate assembly.\n")
            out.write("- Fontes com sinais nesta execução:\n")
            for src, n in sorted(src_counts.items()):
                out.write(f"  - `{src}`: {n} registro(s)\n")

        out.write("\n## 5. Escopo diagnóstico (sem alteração automática de GenBank)\n\n")
        out.write("Este relatório é **diagnostic-only**: ele não altera automaticamente o GenBank e não substitui curadoria manual.\n")

    return md_path, tsv_path
