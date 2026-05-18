from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from .io import infer_format
from .utils import ensure_dir, get_genetic_code, safe_get


def _write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _normalize_organism(config: dict, organism: Optional[str]) -> str:
    if organism:
        return organism
    for path in (["organism"], ["project", "organism"], ["sample"], ["sample_name"]):
        value = safe_get(config, path, None)
        if value:
            return str(value)
    return "Unknown organism"


def _parse_input_record(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Input molecule not found: {path}")
    fmt = infer_format(path)
    return SeqIO.read(str(path), fmt), fmt


def _as_intervals(mask: List[bool], seq_len: int) -> List[Tuple[int, int]]:
    if not any(mask):
        return []
    intervals: List[Tuple[int, int]] = []
    in_run = False
    start = 0
    for i, flagged in enumerate(mask):
        if flagged and not in_run:
            start = i
            in_run = True
        elif not flagged and in_run:
            intervals.append((start, i - 1))
            in_run = False
    if in_run:
        intervals.append((start, seq_len - 1))

    if len(intervals) > 1 and intervals[0][0] == 0 and intervals[-1][1] == seq_len - 1:
        first = intervals.pop(0)
        last = intervals.pop(-1)
        intervals.insert(0, (last[0], first[1] + seq_len))
    return intervals


def detect_at_rich_region(
    seq: str,
    window_size: int,
    min_at_pct: float,
    min_region_len: int,
) -> Optional[Dict[str, Any]]:
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return None
    if window_size <= 0:
        window_size = 200
    if min_region_len <= 0:
        min_region_len = window_size
    if window_size > n:
        window_size = n

    circular = seq + seq[: window_size - 1]
    mask = [False] * n

    at_count = sum(1 for b in circular[:window_size] if b in {"A", "T"})
    if (at_count / window_size) * 100 >= min_at_pct:
        mask[0] = True

    for i in range(1, n):
        left = circular[i - 1]
        right = circular[i + window_size - 1]
        if left in {"A", "T"}:
            at_count -= 1
        if right in {"A", "T"}:
            at_count += 1
        if (at_count / window_size) * 100 >= min_at_pct:
            mask[i] = True

    intervals = _as_intervals(mask, n)
    candidates: List[Dict[str, Any]] = []
    for s, e in intervals:
        length = e - s + 1
        if length < min_region_len:
            continue
        region_seq = (seq + seq)[s : e + 1]
        at_pct = (sum(1 for b in region_seq if b in {"A", "T"}) / len(region_seq)) * 100
        start_1 = (s % n) + 1
        end_1 = (e % n) + 1
        wraps = "yes" if e >= n else "no"
        candidates.append(
            {
                "start": start_1,
                "end": end_1,
                "length": len(region_seq),
                "at_percent": round(at_pct, 3),
                "wraps_around_origin": wraps,
                "score": (round(at_pct, 6), len(region_seq)),
            }
        )

    if not candidates:
        return None
    candidates.sort(key=lambda x: x["score"], reverse=True)
    best = candidates[0]
    best.pop("score", None)
    return best


def _feature_inventory(record) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for feat in record.features:
        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        strand = int(feat.location.strand or 0)
        gene = ";".join(feat.qualifiers.get("gene", [])) or "."
        product = ";".join(feat.qualifiers.get("product", [])) or "."
        notes = ";".join(feat.qualifiers.get("note", [])) or "."
        rows.append(
            {
                "feature_type": feat.type,
                "gene": gene,
                "product": product,
                "start": start,
                "end": end,
                "strand": strand,
                "length_nt": end - start + 1,
                "qualifiers_summary": notes,
            }
        )
    return rows


def run_final_molecule_preparation(
    config: dict,
    root: Path,
    outdir: Path,
    input_molecule: Optional[str] = None,
    organism: Optional[str] = None,
    expected_length: Optional[int] = None,
    topology: str = "circular",
    transl_table: Optional[int] = None,
    at_window_size: int = 500,
    at_min_pct: float = 85.0,
    at_min_len: int = 500,
    existing_annotation: Optional[str] = None,
) -> Path:
    ensure_dir(outdir)

    input_path = Path(input_molecule) if input_molecule else None
    if input_path is None:
        input_path = root / "11_recovered_contig_annotation" / "recovered_contig.best.gb"
        if not input_path.exists():
            gb_candidates = list((root / "11_recovered_contig_annotation").glob("*.gb*"))
            if gb_candidates:
                input_path = gb_candidates[0]
    if not input_path.exists() and existing_annotation:
        input_path = Path(existing_annotation)

    record, in_fmt = _parse_input_record(input_path)
    org = _normalize_organism(config, organism)
    tt = int(transl_table if transl_table is not None else get_genetic_code(config, default=5))

    record.annotations["topology"] = topology
    record.annotations["molecule_type"] = "DNA"
    record.annotations["data_file_division"] = "INV"
    record.annotations.setdefault("keywords", ["."])

    source_feature = SeqFeature(
        location=FeatureLocation(0, len(record.seq), strand=1),
        type="source",
        qualifiers={
            "organism": [org],
            "organelle": ["mitochondrion"],
            "mol_type": ["genomic DNA"],
        },
    )

    non_source_features = [f for f in record.features if f.type != "source"]
    for feat in non_source_features:
        if feat.type == "CDS":
            feat.qualifiers["transl_table"] = [str(tt)]

    at_region = detect_at_rich_region(str(record.seq), at_window_size, at_min_pct, at_min_len)
    if at_region:
        s0 = at_region["start"] - 1
        e0 = at_region["end"]
        if at_region["wraps_around_origin"] == "yes":
            loc = FeatureLocation(s0, len(record.seq), strand=1) + FeatureLocation(0, e0, strand=1)
        else:
            loc = FeatureLocation(s0, e0, strand=1)
        at_feature = SeqFeature(
            location=loc,
            type="misc_feature",
            qualifiers={
                "note": [
                    "A+T-rich control region",
                    "putative mitochondrial control region; not automatically classified as D-loop",
                ]
            },
        )
        non_source_features.append(at_feature)

    record.features = [source_feature] + non_source_features

    fasta_out = outdir / "final_molecule.fasta"
    gb_out = outdir / "final_molecule.gb"
    inv_out = outdir / "final_molecule_feature_inventory.tsv"
    at_out = outdir / "at_rich_region.tsv"
    report_out = outdir / "final_molecule_report.md"

    SeqIO.write(record, str(fasta_out), "fasta")
    SeqIO.write(record, str(gb_out), "genbank")

    inventory = _feature_inventory(record)
    _write_tsv(
        inv_out,
        inventory,
        ["feature_type", "gene", "product", "start", "end", "strand", "length_nt", "qualifiers_summary"],
    )

    at_rows: List[Dict[str, Any]] = []
    if at_region:
        at_rows.append(at_region)
    else:
        at_rows.append({"start": ".", "end": ".", "length": 0, "at_percent": ".", "wraps_around_origin": "no", "status": "not detected"})
    _write_tsv(at_out, at_rows, list(at_rows[0].keys()))

    cds_count = sum(1 for f in record.features if f.type == "CDS")
    rrna_count = sum(1 for f in record.features if f.type == "rRNA")
    trna_count = sum(1 for f in record.features if f.type == "tRNA")

    lines = [
        "# MitoCurator final molecule preparation report\n",
        f"- Input molecule: `{input_path}` ({in_fmt})",
        f"- Organism: {org}",
        f"- Final molecule length: {len(record.seq)} bp",
        f"- Topology set in GenBank LOCUS metadata: {topology}",
        "- Molecule type: genomic DNA (source /mol_type)",
        f"- CDS count: {cds_count}",
        f"- rRNA count: {rrna_count}",
        f"- tRNA count: {trna_count}",
        f"- Total features: {len(record.features)}",
    ]
    if expected_length is not None:
        match = "MATCH" if len(record.seq) == int(expected_length) else "MISMATCH"
        lines.append(f"- Expected length: {expected_length} bp -> {match}")
    if at_region:
        lines.extend([
            f"- A+T-rich/control region: {at_region['start']}..{at_region['end']} (length={at_region['length']} bp, AT={at_region['at_percent']}%)",
            f"- Wrap-around region: {at_region['wraps_around_origin']}",
        ])
    else:
        lines.append("- A+T-rich/control region: not detected with configured thresholds")

    lines.extend([
        "",
        "## Interpretation guardrails\n",
        "- This step prepares a candidate final mitochondrial molecule with explicit circular mitochondrial genomic DNA annotation.",
        "- This step does not by itself prove true biological circularity and does not rule out CNVs.",
        "- Coverage, circular junction, and CNV validation are handled in the next planned module `coverage_cnv_assessment.py`.",
    ])

    report_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return outdir
