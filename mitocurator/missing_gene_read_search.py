from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List

from .utils import ensure_dir, safe_get
from .read_support import resolve_read_sets


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


def _load_candidate_regions(root: Path) -> List[Dict[str, Any]]:
    rows = _read_tsv(root / "08_missing_gene_discovery" / "missing_gene_candidates_pre_mapping.tsv")
    out: List[Dict[str, Any]] = []

    for row in rows:
        if row.get("status") != "candidate_retained":
            continue
        try:
            start = int(row["start"])
            end = int(row["end"])
        except Exception:
            continue

        out.append({
            "gene": row.get("gene", "."),
            "candidate_id": row.get("candidate_id", "."),
            "candidate_rank": row.get("rank", "."),
            "candidate_score": row.get("score", "."),
            "seqid": row.get("seqid", "."),
            "start": start,
            "end": end,
            "strand": row.get("strand", "."),
            "length_nt": max(0, end - start + 1),
            "candidate_length_aa": row.get("length_aa", "."),
            "reference_aa_length": row.get("reference_aa_length", "."),
            "pre_mapping_percent_identity": row.get("percent_identity", "."),
            "pre_mapping_aligned_coverage_reference": row.get("aligned_coverage_reference", "."),
            "pre_mapping_length_ratio": row.get("length_ratio", "."),
        })

    return out


def _summarize_bam_region(bam_path: Path, seqid: str, start1: int, end1: int, min_mapq: int) -> Dict[str, Any]:
    import pysam

    if not bam_path.exists() or not Path(f"{bam_path}.bai").exists():
        return {
            "status": "missing_bam",
            "reads_overlapping": 0,
            "primary_reads": 0,
            "reads_passing_mapq": 0,
            "bases_covered": 0,
            "pct_bases_covered": 0.0,
            "mean_depth": 0.0,
            "min_depth": 0,
            "max_depth": 0,
        }

    start0 = max(0, start1 - 1)
    end0 = max(start0, end1)
    length = max(1, end0 - start0)

    try:
        aln = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as exc:
        return {
            "status": f"bam_open_error:{exc}",
            "reads_overlapping": 0,
            "primary_reads": 0,
            "reads_passing_mapq": 0,
            "bases_covered": 0,
            "pct_bases_covered": 0.0,
            "mean_depth": 0.0,
            "min_depth": 0,
            "max_depth": 0,
        }

    reads_overlapping = 0
    primary_reads = 0
    reads_passing_mapq = 0

    try:
        for read in aln.fetch(seqid, start0, end0):
            reads_overlapping += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            primary_reads += 1
            if int(read.mapping_quality) >= min_mapq:
                reads_passing_mapq += 1
    except ValueError:
        aln.close()
        return {
            "status": "seqid_not_found_in_bam",
            "reads_overlapping": 0,
            "primary_reads": 0,
            "reads_passing_mapq": 0,
            "bases_covered": 0,
            "pct_bases_covered": 0.0,
            "mean_depth": 0.0,
            "min_depth": 0,
            "max_depth": 0,
        }

    depths = [0] * length
    try:
        for col in aln.pileup(seqid, start0, end0, truncate=True, min_mapping_quality=min_mapq):
            pos = col.reference_pos
            if start0 <= pos < end0:
                depths[pos - start0] = col.nsegments
    finally:
        aln.close()

    bases_covered = sum(1 for d in depths if d > 0)
    mean_depth = sum(depths) / length if length else 0.0
    min_depth = min(depths) if depths else 0
    max_depth = max(depths) if depths else 0

    status = "supported" if bases_covered > 0 and reads_passing_mapq > 0 else "no_read_support"

    return {
        "status": status,
        "reads_overlapping": reads_overlapping,
        "primary_reads": primary_reads,
        "reads_passing_mapq": reads_passing_mapq,
        "bases_covered": bases_covered,
        "pct_bases_covered": round(100 * bases_covered / length, 3),
        "mean_depth": round(mean_depth, 3),
        "min_depth": min_depth,
        "max_depth": max_depth,
    }


def _candidate_completeness_class(length_ratio: Any) -> str:
    try:
        ratio = float(length_ratio)
    except Exception:
        return "UNKNOWN_COMPLETENESS"

    if ratio >= 0.90:
        return "NEAR_FULL_LENGTH"
    if ratio >= 0.60:
        return "PARTIAL_HIGH"
    if ratio >= 0.40:
        return "PARTIAL_MEDIUM"
    return "PARTIAL_LOW"


def _combined_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[tuple[str, str], List[Dict[str, Any]]] = {}

    for row in rows:
        grouped.setdefault((row["gene"], row["candidate_id"]), []).append(row)

    out: List[Dict[str, Any]] = []

    for (gene, candidate_id), sub in grouped.items():
        supported = [r for r in sub if r["support_status"] == "supported"]
        n_readsets_evaluated = len(sub)
        n_readsets_supporting = len(supported)

        best_mean_depth = max(float(r["mean_depth"]) for r in sub) if sub else 0.0
        best_pct_cov = max(float(r["pct_bases_covered"]) for r in sub) if sub else 0.0
        total_reads = sum(int(r["reads_passing_mapq"]) for r in sub)

        first = sub[0]
        completeness = _candidate_completeness_class(first.get("pre_mapping_length_ratio", "."))

        if n_readsets_supporting > 0 and best_pct_cov >= 80:
            if completeness == "NEAR_FULL_LENGTH":
                recommendation = "REGION_READ_COVERED_NEAR_FULL_LENGTH_CANDIDATE"
            else:
                recommendation = "REGION_READ_COVERED_PARTIAL_CANDIDATE"
        elif n_readsets_supporting > 0:
            recommendation = "REGION_PARTIALLY_READ_COVERED_PARTIAL_CANDIDATE"
        else:
            recommendation = "REGION_NOT_READ_COVERED"

        final_candidate_decision = "DEFER_TO_TARGETED_CONSENSUS_OR_PROTEIN_READ_SEARCH"

        out.append({
            "gene": gene,
            "candidate_id": candidate_id,
            "candidate_rank": first["candidate_rank"],
            "seqid": first["seqid"],
            "start": first["start"],
            "end": first["end"],
            "strand": first["strand"],
            "length_nt": first["length_nt"],
            "candidate_length_aa": first.get("candidate_length_aa", "."),
            "reference_aa_length": first.get("reference_aa_length", "."),
            "candidate_completeness_class": completeness,
            "n_readsets_evaluated": n_readsets_evaluated,
            "n_readsets_supporting": n_readsets_supporting,
            "readsets_supporting": ",".join(r["read_set"] for r in supported) if supported else ".",
            "total_reads_passing_mapq": total_reads,
            "best_pct_bases_covered": round(best_pct_cov, 3),
            "best_mean_depth": round(best_mean_depth, 3),
            "pre_mapping_percent_identity": first["pre_mapping_percent_identity"],
            "pre_mapping_aligned_coverage_reference": first["pre_mapping_aligned_coverage_reference"],
            "pre_mapping_length_ratio": first["pre_mapping_length_ratio"],
            "recommendation": recommendation,
            "final_candidate_decision": final_candidate_decision,
        })

    # Keep the pre-read candidate rank visible. Read coverage alone cannot
    # select the final missing-gene candidate because all mitochondrial regions
    # may be deeply covered. Final selection is deferred to consensus/protein
    # read-backed evidence.
    out.sort(
        key=lambda r: (
            str(r["gene"]),
            int(float(r["candidate_rank"])) if str(r["candidate_rank"]).replace(".", "", 1).isdigit() else 999999,
        )
    )

    return out


def _write_report(path: Path, combined: List[Dict[str, Any]]) -> None:
    lines = []
    lines.append("# MitoCurator missing-gene read search report\n")
    lines.append("This report summarizes read coverage over provisional missing-gene candidate regions. It does not by itself select the final missing gene annotation.\n")

    if not combined:
        lines.append("No candidate read support records were generated.\n")
    else:
        lines.append("## Ranked candidates\n")
        for row in combined:
            lines.append(
                f"- `{row['gene']}` `{row['candidate_id']}`: "
                f"{row['recommendation']}; "
                f"readsets_supporting={row['n_readsets_supporting']}/{row['n_readsets_evaluated']}; "
                f"length={row.get('candidate_length_aa', '.')}/{row.get('reference_aa_length', '.')} aa; "
                f"class={row.get('candidate_completeness_class', '.')}; "
                f"covered={row['best_pct_bases_covered']}%; "
                f"mean_depth={row['best_mean_depth']}; "
                f"reads={row['total_reads_passing_mapq']}; "
                f"decision={row.get('final_candidate_decision', '.')}"
            )
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def run_missing_gene_read_search(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "09_missing_gene_read_search")

    min_mapq = int(safe_get(config, ["missing_gene_read_search", "min_mapping_quality"], 20))

    candidates = _load_candidate_regions(root)
    read_sets = resolve_read_sets(config)

    per_readset_rows: List[Dict[str, Any]] = []

    for candidate in candidates:
        for rs in read_sets:
            bam = root / "08_read_mapping" / f"{rs['name']}.sorted.bam"
            summary = _summarize_bam_region(
                bam,
                str(candidate["seqid"]),
                int(candidate["start"]),
                int(candidate["end"]),
                min_mapq,
            )

            per_readset_rows.append({
                "gene": candidate["gene"],
                "candidate_id": candidate["candidate_id"],
                "candidate_rank": candidate["candidate_rank"],
                "candidate_score": candidate["candidate_score"],
                "seqid": candidate["seqid"],
                "start": candidate["start"],
                "end": candidate["end"],
                "strand": candidate["strand"],
                "length_nt": candidate["length_nt"],
                "candidate_length_aa": candidate.get("candidate_length_aa", "."),
                "reference_aa_length": candidate.get("reference_aa_length", "."),
                "read_set": rs["name"],
                "read_type": rs["type"],
                "bam": str(bam),
                "support_status": summary["status"],
                "reads_overlapping": summary["reads_overlapping"],
                "primary_reads": summary["primary_reads"],
                "reads_passing_mapq": summary["reads_passing_mapq"],
                "bases_covered": summary["bases_covered"],
                "pct_bases_covered": summary["pct_bases_covered"],
                "mean_depth": summary["mean_depth"],
                "min_depth": summary["min_depth"],
                "max_depth": summary["max_depth"],
                "pre_mapping_percent_identity": candidate["pre_mapping_percent_identity"],
                "pre_mapping_aligned_coverage_reference": candidate["pre_mapping_aligned_coverage_reference"],
                "pre_mapping_length_ratio": candidate["pre_mapping_length_ratio"],
            })

    per_fields = [
        "gene",
        "candidate_id",
        "candidate_rank",
        "candidate_score",
        "seqid",
        "start",
        "end",
        "strand",
        "length_nt",
        "candidate_length_aa",
        "reference_aa_length",
        "read_set",
        "read_type",
        "bam",
        "support_status",
        "reads_overlapping",
        "primary_reads",
        "reads_passing_mapq",
        "bases_covered",
        "pct_bases_covered",
        "mean_depth",
        "min_depth",
        "max_depth",
        "pre_mapping_percent_identity",
        "pre_mapping_aligned_coverage_reference",
        "pre_mapping_length_ratio",
    ]

    combined = _combined_rows(per_readset_rows)

    combined_fields = [
        "gene",
        "candidate_id",
        "candidate_rank",
        "seqid",
        "start",
        "end",
        "strand",
        "length_nt",
        "candidate_length_aa",
        "reference_aa_length",
        "candidate_completeness_class",
        "n_readsets_evaluated",
        "n_readsets_supporting",
        "readsets_supporting",
        "total_reads_passing_mapq",
        "best_pct_bases_covered",
        "best_mean_depth",
        "pre_mapping_percent_identity",
        "pre_mapping_aligned_coverage_reference",
        "pre_mapping_length_ratio",
        "recommendation",
        "final_candidate_decision",
    ]

    _write_tsv(outdir / "missing_gene_candidate_read_support.tsv", per_readset_rows, per_fields)
    _write_tsv(outdir / "missing_gene_candidate_read_support_summary.tsv", combined, combined_fields)
    _write_report(outdir / "missing_gene_read_search_report.md", combined)

    return outdir
