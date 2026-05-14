from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Any

from Bio import SeqIO

from .io import infer_format
from .utils import ensure_dir, safe_get


def _run_shell(cmd: str, log_path: Path) -> None:
    ensure_dir(log_path.parent)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(cmd + "\n\n")
        proc = subprocess.run(
            cmd,
            shell=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {proc.returncode}. See log: {log_path}")


def _resolve_source(config: dict, source: str) -> Any:
    current: Any = config
    for part in source.split("."):
        if isinstance(current, list):
            current = current[int(part)]
        elif isinstance(current, dict):
            current = current[part]
        else:
            raise KeyError(f"Cannot resolve source '{source}' at '{part}'")
    return current


def _as_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(x) for x in value]
    return [str(value)]


def _discover_readsets(config: dict) -> List[Dict[str, Any]]:
    """Discover read sets from read_mapping.read_sets, read_support.read_sets, or reads.

    Supported config patterns:
    1. read_mapping.read_sets:
       - name, type, source
    2. read_support.read_sets:
       - name, type, source
    3. reads as dict:
       reads.hifi: [reads...]
       reads.illumina: {r1: ..., r2: ...}
    4. reads as list of dicts in config.example.yaml style.
    """
    readsets = safe_get(config, ["read_mapping", "read_sets"], None)
    if readsets is None:
        readsets = safe_get(config, ["read_support", "read_sets"], None)

    out: List[Dict[str, Any]] = []

    if isinstance(readsets, list):
        for rs in readsets:
            name = str(rs.get("name"))
            rtype = str(rs.get("type", "unknown"))
            source = rs.get("source")
            resolved = _resolve_source(config, source) if source else rs

            if isinstance(resolved, dict) and {"r1", "r2"}.issubset(resolved.keys()):
                out.append({
                    "name": name,
                    "type": rtype,
                    "r1": str(resolved["r1"]),
                    "r2": str(resolved["r2"]),
                    "reads": [],
                })
            elif isinstance(resolved, dict) and {"read1", "read2"}.issubset(resolved.keys()):
                out.append({
                    "name": name,
                    "type": rtype,
                    "r1": str(resolved["read1"]),
                    "r2": str(resolved["read2"]),
                    "reads": [],
                })
            elif isinstance(resolved, dict) and "reads" in resolved:
                out.append({
                    "name": name,
                    "type": rtype,
                    "r1": "",
                    "r2": "",
                    "reads": _as_list(resolved["reads"]),
                })
            else:
                out.append({
                    "name": name,
                    "type": rtype,
                    "r1": "",
                    "r2": "",
                    "reads": _as_list(resolved),
                })

        return out

    reads = config.get("reads", {})
    if isinstance(reads, dict):
        if "hifi" in reads:
            out.append({
                "name": "hifi",
                "type": "pacbio_hifi",
                "r1": "",
                "r2": "",
                "reads": _as_list(reads.get("hifi")),
            })
        if "ont" in reads:
            out.append({
                "name": "ont",
                "type": "ont",
                "r1": "",
                "r2": "",
                "reads": _as_list(reads.get("ont")),
            })
        if "pacbio_clr" in reads:
            out.append({
                "name": "pacbio_clr",
                "type": "pacbio_clr",
                "r1": "",
                "r2": "",
                "reads": _as_list(reads.get("pacbio_clr")),
            })
        if "illumina" in reads and isinstance(reads["illumina"], dict):
            illumina = reads["illumina"]
            out.append({
                "name": "illumina",
                "type": "illumina_pe",
                "r1": str(illumina.get("r1") or illumina.get("read1") or ""),
                "r2": str(illumina.get("r2") or illumina.get("read2") or ""),
                "reads": [],
            })

    elif isinstance(reads, list):
        for i, rs in enumerate(reads):
            out.append({
                "name": str(rs.get("name", f"readset_{i+1}")),
                "type": str(rs.get("type", "unknown")),
                "r1": str(rs.get("r1") or rs.get("read1") or ""),
                "r2": str(rs.get("r2") or rs.get("read2") or ""),
                "reads": _as_list(rs.get("reads")),
            })

    return [r for r in out if r.get("reads") or (r.get("r1") and r.get("r2"))]


def _minimap_preset(read_type: str) -> str:
    rt = read_type.lower()
    if rt in {"illumina_pe", "illumina_se", "short", "short_reads"}:
        return "sr"
    if rt in {"pacbio_hifi", "hifi", "ccs"}:
        return "map-hifi"
    if rt in {"ont", "nanopore"}:
        return "map-ont"
    if rt in {"pacbio_clr", "clr", "pacbio"}:
        return "map-pb"
    return "map-ont"


def _reference_input(config: dict, root: Path) -> Path:
    preferred = root / "05_refinement" / "refined.gb"
    if preferred.exists():
        return preferred
    return Path(str(safe_get(config, ["input", "mitogenome"])))


def _write_reference_fasta(reference_input: Path, out_fasta: Path) -> Tuple[Path, Path | None]:
    ensure_dir(out_fasta.parent)
    fmt = infer_format(reference_input)

    if fmt == "genbank":
        record = SeqIO.read(str(reference_input), "genbank")
        SeqIO.write(record, str(out_fasta), "fasta")
        return out_fasta, reference_input

    if fmt == "fasta":
        shutil.copyfile(reference_input, out_fasta)
        return out_fasta, None

    raise ValueError(f"Unsupported reference format for read mapping: {reference_input}")


def _feature_name(feature) -> str:
    q = feature.qualifiers
    for key in ["gene", "product", "locus_tag", "label"]:
        if key in q and q[key]:
            return str(q[key][0])
    return feature.type


def _feature_positions(feature) -> List[int]:
    positions: List[int] = []
    parts = getattr(feature.location, "parts", [feature.location])
    for part in parts:
        start = int(part.start)
        end = int(part.end)
        positions.extend(range(start, end))
    return positions


def _read_features(genbank_path: Path | None) -> List[Dict[str, Any]]:
    if genbank_path is None or not genbank_path.exists():
        return []

    record = SeqIO.read(str(genbank_path), "genbank")
    features: List[Dict[str, Any]] = []

    for feature in record.features:
        if feature.type not in {"CDS", "rRNA", "tRNA"}:
            continue
        positions = _feature_positions(feature)
        if not positions:
            continue
        features.append({
            "gene": _feature_name(feature),
            "type": feature.type,
            "start": min(positions) + 1,
            "end": max(positions) + 1,
            "length": len(positions),
            "positions": positions,
        })

    return features


def _write_readsets_tsv(readsets: List[Dict[str, Any]], out_tsv: Path) -> None:
    fieldnames = ["name", "type", "preset", "r1", "r2", "reads"]
    with open(out_tsv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rs in readsets:
            row = {
                "name": rs["name"],
                "type": rs["type"],
                "preset": _minimap_preset(rs["type"]),
                "r1": rs.get("r1", ""),
                "r2": rs.get("r2", ""),
                "reads": ";".join(rs.get("reads", [])),
            }
            writer.writerow(row)


def _idxstats(samtools: str, bam: Path, log_path: Path) -> Tuple[int, int]:
    cmd = f"{samtools} idxstats {bam}"
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"samtools idxstats failed for {bam}. See log: {log_path}")

    mapped = 0
    unmapped = 0
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 4:
            mapped += int(parts[2])
            unmapped += int(parts[3])
    return mapped, unmapped


def _run_depth(samtools: str, bam: Path, out_tsv: Path, readset: str, log_path: Path) -> Dict[str, List[int]]:
    cmd = f"{samtools} depth -aa {bam}"
    result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    log_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"samtools depth failed for {bam}. See log: {log_path}")

    depths_by_seq: Dict[str, List[int]] = {}

    mode = "a" if out_tsv.exists() else "w"
    with open(out_tsv, mode, encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        if mode == "w":
            writer.writerow(["readset", "seqid", "position", "depth"])

        for line in result.stdout.splitlines():
            seqid, pos, depth = line.split("\t")[:3]
            depth_i = int(depth)
            writer.writerow([readset, seqid, pos, depth_i])
            depths_by_seq.setdefault(seqid, []).append(depth_i)

    return depths_by_seq


def _depth_stats(depths_by_seq: Dict[str, List[int]]) -> Dict[str, float]:
    depths = [d for values in depths_by_seq.values() for d in values]
    if not depths:
        return {
            "reference_bases": 0,
            "mean_depth": 0.0,
            "min_depth": 0,
            "max_depth": 0,
            "bases_covered": 0,
            "pct_bases_covered": 0.0,
        }

    ref_bases = len(depths)
    bases_covered = sum(1 for d in depths if d > 0)
    return {
        "reference_bases": ref_bases,
        "mean_depth": sum(depths) / ref_bases,
        "min_depth": min(depths),
        "max_depth": max(depths),
        "bases_covered": bases_covered,
        "pct_bases_covered": 100.0 * bases_covered / ref_bases,
    }


def _write_gene_coverage(
    out_tsv: Path,
    readset: str,
    features: List[Dict[str, Any]],
    depths_by_seq: Dict[str, List[int]],
) -> None:
    mode = "a" if out_tsv.exists() else "w"

    with open(out_tsv, mode, encoding="utf-8", newline="") as fh:
        fieldnames = [
            "readset",
            "gene",
            "type",
            "start",
            "end",
            "length",
            "mean_depth",
            "min_depth",
            "max_depth",
            "bases_covered",
            "pct_bases_covered",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        if mode == "w":
            writer.writeheader()

        if not features:
            return

        # Mitochondrial reference is normally a single sequence; use the first depth vector.
        if not depths_by_seq:
            return
        depth_vector = next(iter(depths_by_seq.values()))

        for feature in features:
            values = [
                depth_vector[pos]
                for pos in feature["positions"]
                if 0 <= pos < len(depth_vector)
            ]
            if values:
                bases_covered = sum(1 for d in values if d > 0)
                mean_depth = sum(values) / len(values)
                min_depth = min(values)
                max_depth = max(values)
                pct = 100.0 * bases_covered / len(values)
            else:
                bases_covered = 0
                mean_depth = 0.0
                min_depth = 0
                max_depth = 0
                pct = 0.0

            writer.writerow({
                "readset": readset,
                "gene": feature["gene"],
                "type": feature["type"],
                "start": feature["start"],
                "end": feature["end"],
                "length": feature["length"],
                "mean_depth": f"{mean_depth:.3f}",
                "min_depth": min_depth,
                "max_depth": max_depth,
                "bases_covered": bases_covered,
                "pct_bases_covered": f"{pct:.3f}",
            })


def run_read_mapping(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "08_read_mapping")
    logs = ensure_dir(outdir / "logs")

    minimap2 = str(safe_get(config, ["tools", "minimap2"], "minimap2"))
    samtools = str(safe_get(config, ["tools", "samtools"], "samtools"))

    minimap_threads = int(safe_get(config, ["read_mapping", "minimap2_threads"], safe_get(config, ["threads", "minimap2"], 8)))
    view_threads = int(safe_get(config, ["read_mapping", "samtools_view_threads"], safe_get(config, ["threads", "samtools_view"], 4)))
    sort_threads = int(safe_get(config, ["read_mapping", "samtools_sort_threads"], safe_get(config, ["threads", "samtools_sort"], 4)))
    min_mapq = int(safe_get(config, ["read_mapping", "min_mapping_quality"], safe_get(config, ["read_support", "min_mapping_quality"], 0)))

    ref_input = _reference_input(config, root)
    ref_fasta, ref_genbank = _write_reference_fasta(ref_input, outdir / "mitogenome_reference.fasta")
    features = _read_features(ref_genbank)

    readsets = _discover_readsets(config)
    if not readsets:
        raise ValueError("No read sets found. Configure reads or read_support.read_sets/read_mapping.read_sets.")

    _write_readsets_tsv(readsets, outdir / "readsets.tsv")

    coverage_by_position = outdir / "coverage_by_position.tsv"
    coverage_by_gene = outdir / "coverage_by_gene.tsv"
    for p in [coverage_by_position, coverage_by_gene]:
        if p.exists():
            p.unlink()

    summary_rows: List[Dict[str, str]] = []

    for rs in readsets:
        name = rs["name"]
        rtype = rs["type"]
        preset = _minimap_preset(rtype)
        bam = outdir / f"{name}.sorted.bam"

        if rtype.lower() == "illumina_pe":
            if not rs.get("r1") or not rs.get("r2"):
                raise ValueError(f"Illumina PE readset '{name}' requires r1 and r2.")
            reads_part = f"{rs['r1']} {rs['r2']}"
        else:
            reads = rs.get("reads", [])
            if not reads:
                raise ValueError(f"Readset '{name}' requires reads.")
            reads_part = " ".join(reads)

        cmd = (
            f"{minimap2} -t {minimap_threads} -ax {preset} {ref_fasta} {reads_part} "
            f"| {samtools} view -@ {view_threads} -b -q {min_mapq} - "
            f"| {samtools} sort -@ {sort_threads} -o {bam} -"
        )

        _run_shell(cmd, logs / f"{name}.mapping.log")
        _run_shell(f"{samtools} index {bam}", logs / f"{name}.index.log")

        mapped, unmapped = _idxstats(samtools, bam, logs / f"{name}.idxstats.log")
        depths_by_seq = _run_depth(
            samtools,
            bam,
            coverage_by_position,
            name,
            logs / f"{name}.depth.log",
        )
        _write_gene_coverage(coverage_by_gene, name, features, depths_by_seq)

        stats = _depth_stats(depths_by_seq)

        summary_rows.append({
            "readset": name,
            "type": rtype,
            "preset": preset,
            "bam": str(bam),
            "mapped_reads": str(mapped),
            "unmapped_reads": str(unmapped),
            "reference_bases": str(int(stats["reference_bases"])),
            "mean_depth": f"{stats['mean_depth']:.3f}",
            "min_depth": str(int(stats["min_depth"])),
            "max_depth": str(int(stats["max_depth"])),
            "bases_covered": str(int(stats["bases_covered"])),
            "pct_bases_covered": f"{stats['pct_bases_covered']:.3f}",
        })

    with open(outdir / "mapping_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "readset",
            "type",
            "preset",
            "bam",
            "mapped_reads",
            "unmapped_reads",
            "reference_bases",
            "mean_depth",
            "min_depth",
            "max_depth",
            "bases_covered",
            "pct_bases_covered",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    return outdir
