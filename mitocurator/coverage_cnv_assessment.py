from __future__ import annotations

import csv
import gzip
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import SeqIO

from .io import infer_format
from .utils import ensure_dir

try:
    import pysam  # type: ignore
except Exception:  # pragma: no cover
    pysam = None


def _write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _classify_depth_ratio(ratio: float, dup: float, strong_dup: float, deletion: float) -> str:
    if ratio > strong_dup:
        return "STRONG_POSSIBLE_DUPLICATION"
    if ratio > dup:
        return "POSSIBLE_DUPLICATION"
    if ratio < deletion:
        return "POSSIBLE_DELETION_OR_LOW_COVERAGE"
    return "NO_CNV_SIGNAL"


def _resolve_input_molecule(root: Path, user_input: Optional[str]) -> Path:
    if user_input:
        p = Path(user_input)
        if not p.exists():
            raise FileNotFoundError(f"Input molecule not found: {p}")
        return p
    candidates = [
        root / "12_final_molecule_preparation" / "final_molecule.gb",
        root / "12_final_molecule_preparation" / "final_molecule.fasta",
        root / "15_applied_curation" / "curated_mitogenome.gb",
        root / "15_applied_curation" / "curated_mitogenome.fasta",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("No input molecule found in expected locations; pass --input-molecule.")


def _resolve_bams(root: Path, input_bam: Optional[Sequence[str]]) -> List[Path]:
    if input_bam:
        out = [Path(p) for p in input_bam]
        miss = [str(p) for p in out if not p.exists()]
        if miss:
            raise FileNotFoundError(f"Missing BAM(s): {', '.join(miss)}")
        return out
    bams = sorted((root / "08_read_mapping").glob("*.sorted.bam"))
    if not bams:
        raise FileNotFoundError("No BAMs found in <outdir>/08_read_mapping/*.sorted.bam; pass --input-bam.")
    return bams


def _make_readset_names(bams: Sequence[Path], names: Optional[Sequence[str]]) -> List[str]:
    if names:
        if len(names) != len(bams):
            raise ValueError("--read-set-name must be provided once per --input-bam.")
        return list(names)
    return [b.stem.replace(".sorted", "") for b in bams]


def _contiguous_feature_parts(feature) -> Iterable[Tuple[int, int]]:
    parts = getattr(feature.location, "parts", [feature.location])
    for part in parts:
        yield int(part.start), int(part.end)


def run_coverage_cnv_assessment(
    config: dict,
    root: Path,
    outdir: Path,
    input_molecule: Optional[str] = None,
    input_bam: Optional[Sequence[str]] = None,
    read_set_name: Optional[Sequence[str]] = None,
    window_size: int = 100,
    step_size: Optional[int] = None,
    duplication_ratio: float = 1.5,
    strong_duplication_ratio: float = 2.0,
    deletion_ratio: float = 0.5,
    junction_flank: int = 500,
    min_mapq: int = 20,
    threads: Optional[int] = None,
) -> Path:
    del config, threads
    if pysam is None:
        raise RuntimeError("pysam is required for BAM-based coverage assessment. Please install pysam and rerun.")

    ensure_dir(outdir)
    molecule = _resolve_input_molecule(root, input_molecule)
    fmt = infer_format(molecule)
    record = SeqIO.read(str(molecule), fmt)
    seqid = record.id
    seq_len = len(record.seq)
    bams = _resolve_bams(root, input_bam)
    readsets = _make_readset_names(bams, read_set_name)

    step = step_size if step_size is not None else window_size
    if window_size <= 0 or step <= 0:
        raise ValueError("window-size and step-size must be positive integers.")

    cov_gz = outdir / "coverage_by_position.tsv.gz"
    windows_out = outdir / "coverage_windows.tsv"
    features_out = outdir / "coverage_by_feature.tsv"
    cnv_out = outdir / "cnv_candidates.tsv"
    junction_out = outdir / "circular_junction_coverage.tsv"
    report_out = outdir / "coverage_cnv_report.md"

    window_rows: List[Dict[str, Any]] = []
    feature_rows: List[Dict[str, Any]] = []
    cnv_rows: List[Dict[str, Any]] = []
    junction_rows: List[Dict[str, Any]] = []
    readset_medians: Dict[str, float] = {}

    with gzip.open(cov_gz, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["read_set", "seqid", "position_1based", "depth"])

        for bam_path, readset in zip(bams, readsets):
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                arr = bam.count_coverage(contig=seqid, start=0, stop=seq_len, quality_threshold=min_mapq)
                depths = [int(arr[0][i] + arr[1][i] + arr[2][i] + arr[3][i]) for i in range(seq_len)]
                for i, d in enumerate(depths, start=1):
                    w.writerow([readset, seqid, i, d])

                gmed = float(median(depths)) if depths else 0.0
                readset_medians[readset] = gmed

                for start0 in range(0, seq_len, step):
                    end0 = min(start0 + window_size, seq_len)
                    segment = depths[start0:end0]
                    if not segment:
                        continue
                    mean_d = sum(segment) / len(segment)
                    med_d = float(median(segment))
                    ratio = (mean_d / gmed) if gmed > 0 else 0.0
                    interp = _classify_depth_ratio(ratio, duplication_ratio, strong_duplication_ratio, deletion_ratio)
                    row = {
                        "read_set": readset, "seqid": seqid,
                        "window_start": start0 + 1, "window_end": end0, "window_len": len(segment),
                        "mean_depth": round(mean_d, 4), "median_depth": round(med_d, 4),
                        "min_depth": min(segment), "max_depth": max(segment),
                        "genome_median_depth": round(gmed, 4), "depth_ratio": round(ratio, 4),
                        "interpretation": interp,
                    }
                    window_rows.append(row)

                if fmt == "genbank":
                    for feat in record.features:
                        ranges = list(_contiguous_feature_parts(feat))
                        vals: List[int] = []
                        for s, e in ranges:
                            vals.extend(depths[s:e])
                        if not vals:
                            continue
                        start = min(s for s, _ in ranges) + 1
                        end = max(e for _, e in ranges)
                        mean_d = sum(vals) / len(vals)
                        med_d = float(median(vals))
                        ratio = (mean_d / gmed) if gmed > 0 else 0.0
                        feature_rows.append({
                            "read_set": readset,
                            "feature_type": feat.type,
                            "gene": ";".join(feat.qualifiers.get("gene", [])) or ".",
                            "product": ";".join(feat.qualifiers.get("product", [])) or ".",
                            "start": start,
                            "end": end,
                            "strand": int(feat.location.strand or 0),
                            "length": len(vals),
                            "mean_depth": round(mean_d, 4),
                            "median_depth": round(med_d, 4),
                            "min_depth": min(vals),
                            "max_depth": max(vals),
                            "genome_median_depth": round(gmed, 4),
                            "depth_ratio": round(ratio, 4),
                            "interpretation": _classify_depth_ratio(ratio, duplication_ratio, strong_duplication_ratio, deletion_ratio),
                        })

                flank = max(1, min(junction_flank, seq_len))
                start_vals = depths[:flank]
                end_vals = depths[seq_len - flank:]
                start_mean = (sum(start_vals) / len(start_vals)) if start_vals else 0.0
                end_mean = (sum(end_vals) / len(end_vals)) if end_vals else 0.0
                start_ratio = (start_mean / gmed) if gmed > 0 else 0.0
                end_ratio = (end_mean / gmed) if gmed > 0 else 0.0
                spanning = 0
                for aln in bam.fetch(seqid, max(0, seq_len - flank - 1), seq_len):
                    if aln.is_unmapped or aln.mapping_quality < min_mapq:
                        continue
                    if (aln.reference_start <= flank) and (aln.reference_end is not None and aln.reference_end >= seq_len - flank):
                        spanning += 1
                if start_ratio >= deletion_ratio and end_ratio >= deletion_ratio:
                    j_interp = "JUNCTION_COVERAGE_COMPATIBLE"
                elif start_ratio < deletion_ratio and end_ratio >= deletion_ratio:
                    j_interp = "LOW_START_FLANK_COVERAGE"
                elif end_ratio < deletion_ratio and start_ratio >= deletion_ratio:
                    j_interp = "LOW_END_FLANK_COVERAGE"
                elif gmed <= 0:
                    j_interp = "JUNCTION_NOT_ASSESSED"
                else:
                    j_interp = "INSUFFICIENT_EVIDENCE"
                junction_rows.append({
                    "read_set": readset, "seqid": seqid, "molecule_length": seq_len, "junction_flank": flank,
                    "start_flank_mean_depth": round(start_mean, 4), "end_flank_mean_depth": round(end_mean, 4),
                    "genome_median_depth": round(gmed, 4), "start_ratio": round(start_ratio, 4), "end_ratio": round(end_ratio, 4),
                    "junction_coverage_interpretation": j_interp, "spanning_read_count": spanning,
                    "comment": "Spanning-read count from linear reference BAM may underestimate true circular junction support.",
                })

    # CNV from windows then merge adjacent same-signal per readset
    for readset in readsets:
        rows = [r for r in window_rows if r["read_set"] == readset and r["interpretation"] != "NO_CNV_SIGNAL"]
        rows.sort(key=lambda x: int(x["window_start"]))
        cur = None
        for r in rows:
            if cur and r["interpretation"] == cur["signal"] and int(r["window_start"]) <= int(cur["end"]) + 1:
                cur["end"] = int(r["window_end"])
                cur["length"] = cur["end"] - cur["start"] + 1
                cur["_means"].append(float(r["mean_depth"]))
                cur["_meds"].append(float(r["median_depth"]))
                cur["_ratios"].append(float(r["depth_ratio"]))
            else:
                if cur:
                    cnv_rows.append({k: v for k, v in cur.items() if not k.startswith("_")})
                cur = {
                    "read_set": readset, "seqid": seqid, "start": int(r["window_start"]), "end": int(r["window_end"]),
                    "length": int(r["window_len"]), "mean_depth": float(r["mean_depth"]), "median_depth": float(r["median_depth"]),
                    "genome_median_depth": float(r["genome_median_depth"]), "depth_ratio": float(r["depth_ratio"]),
                    "signal": r["interpretation"], "evidence_type": "window_merged", "comment": "Merged adjacent windows with same CNV signal.",
                    "_means": [float(r["mean_depth"])], "_meds": [float(r["median_depth"])], "_ratios": [float(r["depth_ratio"])],
                }
        if cur:
            cur["mean_depth"] = round(sum(cur["_means"]) / len(cur["_means"]), 4)
            cur["median_depth"] = round(sum(cur["_meds"]) / len(cur["_meds"]), 4)
            cur["depth_ratio"] = round(sum(cur["_ratios"]) / len(cur["_ratios"]), 4)
            cnv_rows.append({k: v for k, v in cur.items() if not k.startswith("_")})

    _write_tsv(windows_out, window_rows, ["read_set", "seqid", "window_start", "window_end", "window_len", "mean_depth", "median_depth", "min_depth", "max_depth", "genome_median_depth", "depth_ratio", "interpretation"])

    feature_fields = ["read_set", "feature_type", "gene", "product", "start", "end", "strand", "length", "mean_depth", "median_depth", "min_depth", "max_depth", "genome_median_depth", "depth_ratio", "interpretation"]
    _write_tsv(features_out, feature_rows, feature_fields)
    _write_tsv(cnv_out, cnv_rows, ["read_set", "seqid", "start", "end", "length", "mean_depth", "median_depth", "genome_median_depth", "depth_ratio", "signal", "evidence_type", "comment"])
    _write_tsv(junction_out, junction_rows, ["read_set", "seqid", "molecule_length", "junction_flank", "start_flank_mean_depth", "end_flank_mean_depth", "genome_median_depth", "start_ratio", "end_ratio", "junction_coverage_interpretation", "spanning_read_count", "comment"])

    feature_note = "Feature-level coverage computed from GenBank features." if fmt == "genbank" else "Feature-level coverage requires GenBank input; FASTA input produced header-only feature table."
    lines = [
        "# Coverage and CNV assessment report\n",
        f"- Input molecule: `{molecule}`",
        f"- Molecule length: {seq_len} bp",
        f"- Input format: {fmt}",
        f"- BAM/read sets evaluated: {', '.join(f'{rs}:{bp}' for rs, bp in zip(readsets, bams))}",
        "- Genome-wide median depth per read set:",
    ]
    for rs in readsets:
        lines.append(f"  - {rs}: {readset_medians.get(rs, 0.0):.4f}")
    lines.extend([
        f"- Window CNV signal count (non-NO_CNV_SIGNAL): {sum(1 for r in window_rows if r['interpretation'] != 'NO_CNV_SIGNAL')}",
        f"- Feature CNV signal count (non-NO_CNV_SIGNAL): {sum(1 for r in feature_rows if r['interpretation'] != 'NO_CNV_SIGNAL')}",
        "- Circular junction assessment summary:",
    ])
    for r in junction_rows:
        lines.append(f"  - {r['read_set']}: {r['junction_coverage_interpretation']} (start_ratio={r['start_ratio']}, end_ratio={r['end_ratio']}, spanning_read_count={r['spanning_read_count']})")
    lines.extend([
        f"- Feature table note: {feature_note}",
        "",
        "## Interpretation guardrails",
        "- This module evaluates coverage compatibility and CNV signals but does not by itself prove circularity.",
        "- Spanning-read evidence should be interpreted cautiously when reads were mapped to a linearized reference.",
        "- Candidate duplications/deletions require manual review and, ideally, orthogonal evidence.",
    ])
    report_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return outdir
