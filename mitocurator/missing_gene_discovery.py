from __future__ import annotations

import csv
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from .utils import ensure_dir, safe_get


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


def _num(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _score_candidate(candidate: Dict[str, str], similarity: Dict[str, str] | None) -> float:
    score = 0.0

    hint = str(candidate.get("decision_hint", "") or "")
    if "STRONG_CANDIDATE_REGION" in hint:
        score += 10000
    if "LOW_PRIORITY" in hint:
        score -= 1000

    if str(candidate.get("internal_stop_count", "999")) == "0":
        score += 5000

    score += _num(candidate.get("length_aa", "0")) * 5

    if similarity:
        score += _num(similarity.get("aligned_coverage_reference", "0")) * 25
        score += _num(similarity.get("percent_identity", "0")) * 10
        score += _num(similarity.get("aligned_coverage_candidate", "0")) * 2
        score += _num(similarity.get("length_ratio", "0")) * 1000

    return score


def _missing_cds_genes(root: Path) -> List[str]:
    inventory = _read_tsv(root / "06_annotation_assessment" / "annotation_gene_inventory.tsv")
    return sorted({
        row.get("gene", "")
        for row in inventory
        if row.get("type") == "CDS"
        and row.get("annotation_status") == "MISSING_OR_INCOMPLETE"
    })


def _discover_missing_gene_candidates(root: Path, config: dict) -> List[Dict[str, Any]]:
    candidates = _read_tsv(root / "05_refinement" / "missing_gene_candidates.tsv")
    similarities = _read_tsv(root / "05_refinement" / "reference_similarity_candidates.tsv")

    sim_by_key = {
        (row.get("gene", ""), row.get("candidate_id", "")): row
        for row in similarities
    }

    min_ref_cov = float(safe_get(config, ["missing_gene_discovery", "min_aligned_coverage_reference"], 25.0))
    min_pid = float(safe_get(config, ["missing_gene_discovery", "min_percent_identity"], 25.0))
    max_internal_stops = int(safe_get(config, ["missing_gene_discovery", "max_internal_stop_count"], 0))
    keep_low_priority = bool(safe_get(config, ["missing_gene_discovery", "keep_low_priority"], True))

    rows: List[Dict[str, Any]] = []

    for gene in _missing_cds_genes(root):
        gene_rows = [r for r in candidates if r.get("gene") == gene]

        if not gene_rows:
            rows.append({
                "gene": gene,
                "candidate_id": ".",
                "status": "no_candidate_found",
                "rank": ".",
                "score": ".",
                "seqid": ".",
                "start": ".",
                "end": ".",
                "strand": ".",
                "length_nt": ".",
                "length_aa": ".",
                "internal_stop_count": ".",
                "terminal_stop": ".",
                "candidate_decision_hint": ".",
                "reference_gene": ".",
                "reference_aa_length": ".",
                "length_ratio": ".",
                "percent_identity": ".",
                "aligned_coverage_reference": ".",
                "aligned_coverage_candidate": ".",
                "candidate_class": "missing_without_candidate",
                "comment": "No pre-read-mapping candidate found.",
            })
            continue

        scored = []
        for cand in gene_rows:
            sim = sim_by_key.get((gene, cand.get("candidate_id", "")))

            internal_stops = int(_num(cand.get("internal_stop_count", "999")))
            if internal_stops > max_internal_stops:
                status = "filtered_internal_stop_count"
            elif sim and _num(sim.get("aligned_coverage_reference", "0")) < min_ref_cov:
                status = "filtered_low_reference_coverage"
            elif sim and _num(sim.get("percent_identity", "0")) < min_pid:
                status = "filtered_low_percent_identity"
            elif (not keep_low_priority) and "LOW_PRIORITY" in str(cand.get("decision_hint", "")):
                status = "filtered_low_priority"
            else:
                status = "candidate_retained"

            score = _score_candidate(cand, sim)
            scored.append((score, status, cand, sim))

        scored.sort(key=lambda x: x[0], reverse=True)

        for rank, (score, status, cand, sim) in enumerate(scored, start=1):
            rows.append({
                "gene": gene,
                "candidate_id": cand.get("candidate_id", "."),
                "status": status,
                "rank": rank,
                "score": f"{score:.3f}",
                "seqid": cand.get("seqid", "."),
                "start": cand.get("start", "."),
                "end": cand.get("end", "."),
                "strand": cand.get("strand", "."),
                "length_nt": cand.get("length_nt", "."),
                "length_aa": cand.get("length_aa", "."),
                "internal_stop_count": cand.get("internal_stop_count", "."),
                "terminal_stop": cand.get("terminal_stop", "."),
                "candidate_decision_hint": cand.get("decision_hint", "."),
                "reference_gene": sim.get("reference_gene", ".") if sim else ".",
                "reference_aa_length": sim.get("reference_aa_length", ".") if sim else ".",
                "length_ratio": sim.get("length_ratio", ".") if sim else ".",
                "percent_identity": sim.get("percent_identity", ".") if sim else ".",
                "aligned_coverage_reference": sim.get("aligned_coverage_reference", ".") if sim else ".",
                "aligned_coverage_candidate": sim.get("aligned_coverage_candidate", ".") if sim else ".",
                "candidate_class": "pre_read_mapping_candidate",
                "comment": "Candidate retained for downstream read-backed validation; not applied as final CDS.",
            })

    return rows


def _write_bed(path: Path, rows: List[Dict[str, Any]]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as out:
        for row in rows:
            if row.get("status") != "candidate_retained":
                continue
            try:
                start0 = int(row["start"]) - 1
                end = int(row["end"])
            except Exception:
                continue
            name = f"{row.get('gene', '.')}_{row.get('candidate_id', '.')}"
            out.write(f"{row.get('seqid', '.')}\t{start0}\t{end}\t{name}\t{row.get('score', '.')}\t{row.get('strand', '.')}\n")


def _add_candidate_misc_features(record, rows: List[Dict[str, Any]], genetic_code: int):
    curated = deepcopy(record)

    for row in rows:
        if row.get("status") != "candidate_retained":
            continue

        try:
            start = int(row["start"])
            end = int(row["end"])
        except Exception:
            continue

        strand = -1 if str(row.get("strand", "+")) in {"-", "-1"} else 1
        loc = FeatureLocation(start - 1, end, strand=strand)

        note = (
            "MitoCurator: candidate feature for missing gene discovery. "
            "This is not a final CDS annotation; downstream read evidence is required."
        )

        feature = SeqFeature(
            location=loc,
            type="misc_feature",
            qualifiers={
                "gene": [row.get("gene", ".")],
                "label": [f"{row.get('gene', '.')}_{row.get('candidate_id', '.')}"],
                "note": [
                    note,
                    f"candidate_id={row.get('candidate_id', '.')}",
                    f"rank={row.get('rank', '.')}",
                    f"score={row.get('score', '.')}",
                    f"pid={row.get('percent_identity', '.')}",
                    f"aligned_coverage_reference={row.get('aligned_coverage_reference', '.')}",
                    f"length_ratio={row.get('length_ratio', '.')}",
                    f"transl_table={genetic_code}",
                ],
            },
        )

        curated.features.append(feature)

    curated.features.sort(key=lambda f: int(f.location.start))
    curated.description = curated.description + " | MitoCurator missing-gene candidate discovery"
    curated.annotations = deepcopy(record.annotations)
    curated.annotations.setdefault("molecule_type", "DNA")

    return curated


def _write_report(path: Path, rows: List[Dict[str, Any]]) -> None:
    retained = [r for r in rows if r.get("status") == "candidate_retained"]
    genes = sorted({r.get("gene", ".") for r in rows})

    lines = []
    lines.append("# MitoCurator missing-gene discovery report\n")
    lines.append("This report records candidate regions for missing/incomplete CDSs before read-backed validation.\n")
    lines.append("## Summary\n")
    lines.append(f"- Missing CDS genes evaluated: {len(genes)}")
    lines.append(f"- Candidate records: {len(rows)}")
    lines.append(f"- Retained candidates: {len(retained)}\n")

    lines.append("## Retained candidates\n")
    if retained:
        for row in retained:
            lines.append(
                f"- `{row['gene']}` `{row['candidate_id']}` rank={row['rank']} "
                f"{row['start']}-{row['end']} strand={row['strand']} "
                f"score={row['score']} pid={row['percent_identity']} "
                f"cov_ref={row['aligned_coverage_reference']} ratio={row['length_ratio']}"
            )
    else:
        lines.append("- No candidates retained.")
    lines.append("")

    lines.append("These candidates are provisional. Final gene/CDS application must be based on downstream read-backed evidence.")
    path.write_text("\n".join(lines), encoding="utf-8")


def run_missing_gene_discovery(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "08_missing_gene_discovery")

    input_gb = root / "05_refinement" / "refined.gb"
    if not input_gb.exists():
        raise FileNotFoundError(f"Missing refined GenBank: {input_gb}")

    genetic_code = int(
        safe_get(
            config,
            ["missing_gene_discovery", "genetic_code"],
            safe_get(config, ["mitofinder", "organism_code"], 5),
        )
    )

    record = SeqIO.read(str(input_gb), "genbank")
    rows = _discover_missing_gene_candidates(root, config)
    candidate_record = _add_candidate_misc_features(record, rows, genetic_code)

    fields = [
        "gene",
        "candidate_id",
        "status",
        "rank",
        "score",
        "seqid",
        "start",
        "end",
        "strand",
        "length_nt",
        "length_aa",
        "internal_stop_count",
        "terminal_stop",
        "candidate_decision_hint",
        "reference_gene",
        "reference_aa_length",
        "length_ratio",
        "percent_identity",
        "aligned_coverage_reference",
        "aligned_coverage_candidate",
        "candidate_class",
        "comment",
    ]

    _write_tsv(outdir / "missing_gene_candidates_pre_mapping.tsv", rows, fields)
    _write_bed(outdir / "missing_gene_candidate_regions.bed", rows)
    SeqIO.write(candidate_record, str(outdir / "missing_gene_candidate_features.gb"), "genbank")
    SeqIO.write(candidate_record, str(outdir / "missing_gene_candidate_features.fasta"), "fasta")
    _write_report(outdir / "missing_gene_discovery_report.md", rows)

    return outdir
