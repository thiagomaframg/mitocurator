from __future__ import annotations

import csv
import gzip
import random
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

from Bio import SeqIO

from .io import infer_format
from .utils import ensure_dir, safe_get
from .read_support import resolve_read_sets


LONG_READ_TYPES = {"pacbio_hifi", "hifi", "pacbio_clr", "clr", "ont", "nanopore"}


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


def _as_int(value: Any, default: int) -> int:
    try:
        if isinstance(value, dict):
            return default
        if value is None:
            return default
        return int(value)
    except Exception:
        return default


def _resolve_threads(config: dict) -> int:
    value = safe_get(config, ["missing_gene_recovery", "threads"], None)
    if value is not None:
        return _as_int(value, 8)

    for key in ["missing_gene_recovery", "minimap2", "bwa", "samtools"]:
        value = safe_get(config, ["threads", key], None)
        if value is not None:
            return _as_int(value, 8)

    value = safe_get(config, ["threads"], None)
    return _as_int(value, 8)


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


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _open_write_text(path: Path):
    ensure_dir(path.parent)
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "w", encoding="utf-8")


def _fastq_iter(path: Path):
    with _open_text(path) as fh:
        while True:
            name = fh.readline()
            if not name:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            yield name.rstrip(), seq.rstrip(), plus.rstrip(), qual.rstrip()


def _read_id(header: str) -> str:
    h = header[1:] if header.startswith("@") or header.startswith(">") else header
    return h.split()[0]


def _fastq_stats(paths: List[Path]) -> Tuple[int, int]:
    n = 0
    b = 0
    for path in paths:
        for _, seq, _, _ in _fastq_iter(path):
            n += 1
            b += len(seq)
    return n, b


def _write_reads_by_id(input_paths: List[Path], output_path: Path, keep_ids: Set[str]) -> Tuple[int, int]:
    n = 0
    b = 0
    written: Set[str] = set()

    with _open_write_text(output_path) as out:
        for path in input_paths:
            for name, seq, plus, qual in _fastq_iter(path):
                rid = _read_id(name)
                if rid not in keep_ids or rid in written:
                    continue
                out.write(f"{name}\n{seq}\n{plus}\n{qual}\n")
                written.add(rid)
                n += 1
                b += len(seq)

    return n, b


def _missing_cds_genes(root: Path) -> List[str]:
    inv = _read_tsv(root / "06_annotation_assessment" / "annotation_gene_inventory.tsv")
    return sorted({
        row.get("gene", "")
        for row in inv
        if row.get("type") == "CDS"
        and row.get("annotation_status") == "MISSING_OR_INCOMPLETE"
    })


def _extract_reference_genes(config: dict, genes: List[str], out_fasta: Path) -> Path:
    ref_gb = Path(str(safe_get(config, ["mitofinder", "reference_gb"], "")))
    if not ref_gb.exists():
        raise FileNotFoundError(f"Reference GenBank not found: {ref_gb}")

    records = []
    ref = SeqIO.read(str(ref_gb), "genbank")

    for gene in genes:
        found = False
        for feat in ref.features:
            if feat.type != "CDS":
                continue
            if gene not in feat.qualifiers.get("gene", []):
                continue
            seq = feat.extract(ref.seq)
            rec = SeqIO.SeqRecord(
                seq,
                id=gene,
                name=gene,
                description=f"{gene} reference CDS extracted from {ref_gb.name}",
            )
            records.append(rec)
            found = True
            break
        if not found:
            records.append(
                SeqIO.SeqRecord(
                    ref.seq[:0],
                    id=gene,
                    name=gene,
                    description=f"{gene} not found in {ref_gb.name}",
                )
            )

    ensure_dir(out_fasta.parent)
    SeqIO.write(records, str(out_fasta), "fasta")
    return out_fasta


def _write_mitogenome_reference(root: Path, out_fasta: Path) -> Tuple[Path, int]:
    ref_gb = root / "05_refinement" / "refined.gb"
    if not ref_gb.exists():
        raise FileNotFoundError(f"Missing refined GenBank: {ref_gb}")

    record = SeqIO.read(str(ref_gb), "genbank")
    ensure_dir(out_fasta.parent)
    SeqIO.write(record, str(out_fasta), "fasta")
    return out_fasta, len(record.seq)


def _parse_paf(path: Path) -> List[Dict[str, Any]]:
    rows = []
    if not path.exists():
        return rows

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            qname = p[0]
            qlen = int(p[1])
            qstart = int(p[2])
            qend = int(p[3])
            strand = p[4]
            tname = p[5]
            tlen = int(p[6])
            tstart = int(p[7])
            tend = int(p[8])
            matches = int(p[9])
            aln_len = int(p[10])
            mapq = int(p[11])
            rows.append({
                "read_id": qname,
                "query_length": qlen,
                "query_aligned_length": qend - qstart,
                "target": tname,
                "target_length": tlen,
                "target_aligned_length": tend - tstart,
                "strand": strand,
                "matches": matches,
                "alignment_length": aln_len,
                "mapq": mapq,
                "query_coverage": round(100 * (qend - qstart) / qlen, 3) if qlen else 0,
                "target_coverage": round(100 * (tend - tstart) / tlen, 3) if tlen else 0,
                "identity": round(100 * matches / aln_len, 3) if aln_len else 0,
            })

    return rows






def _select_mitogenome_read_ids_for_target_coverage(
    paf_rows: List[Dict[str, Any]],
    min_mapq: int,
    min_identity: float,
    min_alignment_length: int,
    max_reads: int,
) -> tuple[Set[str], int, int, float, str]:
    """Select top high-confidence mitogenome-like reads.

    This reproduces the M. capixaba HiFi empirical strategy: keep the strongest
    mitogenome-like reads after strict PAF filters, then take the top N reads
    for local recovery assembly. In the M. capixaba case, N=300.
    """
    best_by_read: Dict[str, Dict[str, Any]] = {}

    for row in paf_rows:
        if int(row["mapq"]) < min_mapq:
            continue
        if float(row["identity"]) < min_identity:
            continue
        if int(row["alignment_length"]) < min_alignment_length:
            continue

        rid = str(row["read_id"])
        old = best_by_read.get(rid)
        if old is None or int(row["alignment_length"]) > int(old["alignment_length"]):
            best_by_read[rid] = row

    rows = list(best_by_read.values())
    rows.sort(
        key=lambda r: (
            int(r["alignment_length"]),
            float(r["identity"]),
            int(r["mapq"]),
        ),
        reverse=True,
    )

    selected_rows = rows[:max_reads] if max_reads > 0 else rows
    selected_ids = {str(row["read_id"]) for row in selected_rows}

    mean_passing_read_length = (
        sum(int(r["query_length"]) for r in rows) / len(rows)
        if rows else 0.0
    )

    selection_note = f"top_{max_reads}_of_{len(rows)}_passing_reads" if max_reads > 0 else f"all_{len(rows)}_passing_reads"

    return selected_ids, len(rows), len(selected_rows), mean_passing_read_length, selection_note

def _select_gene_read_ids_from_paf(
    paf_rows: List[Dict[str, Any]],
    min_mapq: int,
    min_identity: float,
    min_target_aligned_fraction: float,
) -> Set[str]:
    """Select reads containing a missing reference gene/CDS.

    The length threshold is dynamic. For the M. capixaba ND2 case, the manual
    cutoff of ~700 nt corresponds to ~70% of the 982 nt ND2 reference. The same
    fraction is used for any missing CDS.
    """
    keep = set()

    for row in paf_rows:
        if int(row["mapq"]) < min_mapq:
            continue
        if float(row["identity"]) < min_identity:
            continue

        target_length = int(row["target_length"])
        min_target_aligned_length = int(round(target_length * min_target_aligned_fraction))

        if int(row["target_aligned_length"]) < min_target_aligned_length:
            continue

        keep.add(str(row["read_id"]))

    return keep

def _flye_read_option(read_type: str) -> str:
    rtype = str(read_type).lower()
    if rtype in {"pacbio_hifi", "hifi"}:
        return "--pacbio-hifi"
    if rtype in {"ont", "nanopore"}:
        return "--nano-raw"
    if rtype in {"pacbio_clr", "clr"}:
        return "--pacbio-raw"
    return "--pacbio-hifi"


def _run_flye_for_pool(
    config: dict,
    read_set_name: str,
    read_type: str,
    pool_fastq: Path,
    selected_dir: Path,
    threads: int,
    genome_size: str,
    log_path: Path,
) -> Path:
    flye = str(safe_get(config, ["tools", "flye"], "flye"))
    read_option = _flye_read_option(read_type)
    outdir = selected_dir / f"{read_set_name}.flye_missing_gene_recovery"

    cmd = (
        f"{flye} {read_option} {pool_fastq} "
        f"--genome-size {genome_size} "
        f"--threads {threads} "
        f"--out-dir {outdir}"
    )

    _run_shell(cmd, log_path)
    return outdir


def _long_readsets(config: dict) -> List[Dict[str, Any]]:
    out = []
    for rs in resolve_read_sets(config):
        rtype = str(rs.get("type", "")).lower()
        if rtype not in LONG_READ_TYPES:
            continue
        reads = rs.get("reads")
        if not reads:
            continue
        if isinstance(reads, str):
            reads = [reads]
        out.append({**rs, "reads": [Path(str(x)) for x in reads]})
    return out


def run_missing_gene_recovery(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "09_missing_gene_recovery")
    logs = ensure_dir(outdir / "logs")
    selected_dir = ensure_dir(outdir / "selected_reads")

    genes = _missing_cds_genes(root)
    ref_gene_fasta = _extract_reference_genes(
        config,
        genes,
        outdir / "missing_reference_genes.fasta",
    )
    mito_fasta, mito_len = _write_mitogenome_reference(
        root,
        outdir / "mitogenome_reference.fasta",
    )

    threads = _resolve_threads(config)
    target_cov = float(safe_get(config, ["missing_gene_recovery", "target_mitogenome_coverage"], 400))
    min_mapq = int(safe_get(config, ["missing_gene_recovery", "min_mapq"], 30))

    # Missing-gene/CDS reads. The default 0.70 reproduces the ND2 manual
    # threshold (~700 nt over a 982 nt reference) without hard-coding ND2.
    gene_min_identity = float(safe_get(config, ["missing_gene_recovery", "gene_min_identity"], 0))
    gene_min_target_aligned_fraction = float(
        safe_get(config, ["missing_gene_recovery", "gene_min_target_aligned_fraction"], 0.713)
    )

    # Mitogenome-like reads used as the backbone local assembly pool.
    mito_min_identity = float(safe_get(config, ["missing_gene_recovery", "mito_min_identity"], 0))
    mito_min_alignment_length = int(
        safe_get(config, ["missing_gene_recovery", "mito_min_alignment_length"], 10000)
    )
    hifi_mitogenome_reads = int(safe_get(config, ["missing_gene_recovery", "hifi_mitogenome_reads"], 300))
    mito_min_reads = int(safe_get(config, ["missing_gene_recovery", "mito_min_reads"], 300))
    seed = int(safe_get(config, ["missing_gene_recovery", "random_seed"], 42))
    rng = random.Random(seed)

    run_flye = bool(safe_get(config, ["missing_gene_recovery", "run_flye"], False))
    flye_genome_size = str(safe_get(config, ["missing_gene_recovery", "flye_genome_size"], "20k"))

    readsets = _long_readsets(config)

    gene_hit_rows = []
    mito_hit_rows = []
    pool_rows = []

    for rs in readsets:
        name = rs["name"]
        rtype = str(rs["type"]).lower()
        preset = "map-hifi" if rtype in {"pacbio_hifi", "hifi"} else ("map-ont" if rtype in {"ont", "nanopore"} else "map-pb")
        read_args = " ".join(str(p) for p in rs["reads"])

        gene_paf = outdir / f"{name}.missing_gene_vs_reads.paf"
        mito_paf = outdir / f"{name}.mitogenome_vs_reads.paf"

        _run_shell(
            f"minimap2 -t {threads} -x {preset} {ref_gene_fasta} {read_args} > {gene_paf}",
            logs / f"{name}.missing_gene_vs_reads.log",
        )
        _run_shell(
            f"minimap2 -t {threads} -x {preset} {mito_fasta} {read_args} > {mito_paf}",
            logs / f"{name}.mitogenome_vs_reads.log",
        )

        gene_paf_rows = _parse_paf(gene_paf)
        mito_paf_rows = _parse_paf(mito_paf)

        for row in gene_paf_rows:
            gene_hit_rows.append({"read_set": name, **row})
        for row in mito_paf_rows:
            mito_hit_rows.append({"read_set": name, **row})

        gene_ids = _select_gene_read_ids_from_paf(
            gene_paf_rows,
            min_mapq=min_mapq,
            min_identity=gene_min_identity,
            min_target_aligned_fraction=gene_min_target_aligned_fraction,
        )
        mito_ids_all, mito_ids_passing_filters, estimated_target_mito_reads, mean_passing_mito_read_len, mito_selection_note = (
            _select_mitogenome_read_ids_for_target_coverage(
                mito_paf_rows,
                min_mapq=min_mapq,
                min_identity=mito_min_identity,
                min_alignment_length=mito_min_alignment_length,
                max_reads=hifi_mitogenome_reads,
            )
        )

        n_reads, n_bases = _fastq_stats(rs["reads"])
        mean_read_len = n_bases / n_reads if n_reads else 0

        mito_ids = sorted(mito_ids_all)
        final_ids = set(mito_ids) | gene_ids

        out_fastq = selected_dir / f"{name}.missing_gene_recovery_pool.fastq.gz"
        written_reads, written_bases = _write_reads_by_id(rs["reads"], out_fastq, final_ids)

        flye_script = selected_dir / f"{name}.missing_gene_recovery.flye.sh"
        flye_outdir = selected_dir / f"{name}.flye_missing_gene_recovery"
        flye_read_option = _flye_read_option(rtype)

        flye_script.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n\n"
            f"{str(safe_get(config, ['tools', 'flye'], 'flye'))} {flye_read_option} {out_fastq} "
            f"--genome-size {flye_genome_size} "
            f"--threads {threads} --out-dir {flye_outdir}\n",
            encoding="utf-8",
        )
        flye_script.chmod(0o755)

        flye_status = "not_run"
        flye_assembly = flye_outdir / "assembly.fasta"

        if run_flye:
            try:
                _run_flye_for_pool(
                    config=config,
                    read_set_name=name,
                    read_type=rtype,
                    pool_fastq=out_fastq,
                    selected_dir=selected_dir,
                    threads=threads,
                    genome_size=flye_genome_size,
                    log_path=logs / f"{name}.flye_missing_gene_recovery.log",
                )
                flye_status = "completed" if flye_assembly.exists() else "completed_no_assembly"
            except Exception as exc:
                flye_status = f"failed:{exc}"

        pool_rows.append({
            "read_set": name,
            "read_type": rtype,
            "input_reads": n_reads,
            "input_bases": n_bases,
            "mitogenome_length": mito_len,
            "target_mitogenome_coverage": target_cov,
            "mean_read_length": round(mean_read_len, 3),
            "mitogenome_read_ids_passing_filters": mito_ids_passing_filters,
            "mitogenome_target_reads_requested": hifi_mitogenome_reads,
            "mitogenome_read_ids_after_subsampling": len(mito_ids),
            "mean_passing_mitogenome_read_length": round(mean_passing_mito_read_len, 3),
            "missing_gene_read_ids": len(gene_ids),
            "gene_min_target_aligned_fraction": gene_min_target_aligned_fraction,
            "mito_min_alignment_length": mito_min_alignment_length,
            "min_mapq": min_mapq,
            "final_pool_read_ids": len(final_ids),
            "written_reads": written_reads,
            "written_bases": written_bases,
            "approx_pool_coverage": round(written_bases / mito_len, 3) if mito_len else 0,
            "output_fastq": str(out_fastq),
            "flye_script": str(flye_script),
            "run_flye": run_flye,
            "flye_genome_size": flye_genome_size,
            "flye_status": flye_status,
            "flye_outdir": str(flye_outdir),
            "flye_assembly": str(flye_assembly) if flye_assembly.exists() else ".",
        })

    _write_tsv(
        outdir / "missing_gene_read_hits.tsv",
        gene_hit_rows,
        [
            "read_set",
            "read_id",
            "query_length",
            "query_aligned_length",
            "target",
            "target_length",
            "target_aligned_length",
            "strand",
            "matches",
            "alignment_length",
            "mapq",
            "query_coverage",
            "target_coverage",
            "identity",
        ],
    )

    _write_tsv(
        outdir / "mitogenome_read_hits.tsv",
        mito_hit_rows,
        [
            "read_set",
            "read_id",
            "query_length",
            "query_aligned_length",
            "target",
            "target_length",
            "target_aligned_length",
            "strand",
            "matches",
            "alignment_length",
            "mapq",
            "query_coverage",
            "target_coverage",
            "identity",
        ],
    )

    _write_tsv(
        outdir / "selected_read_pool.tsv",
        pool_rows,
        [
            "read_set",
            "read_type",
            "input_reads",
            "input_bases",
            "mitogenome_length",
            "target_mitogenome_coverage",
            "mean_read_length",
            "mitogenome_read_ids_passing_filters",
            "mitogenome_target_reads_requested",
            "mitogenome_read_ids_after_subsampling",
            "mean_passing_mitogenome_read_length",
            "missing_gene_read_ids",
            "gene_min_target_aligned_fraction",
            "mito_min_alignment_length",
            "min_mapq",
            "final_pool_read_ids",
            "written_reads",
            "written_bases",
            "approx_pool_coverage",
            "output_fastq",
            "flye_script",
            "run_flye",
            "flye_genome_size",
            "flye_status",
            "flye_outdir",
            "flye_assembly",
        ],
    )

    lines = []
    lines.append("# MitoCurator missing-gene recovery report\n")
    lines.append("This step follows the direct recovery strategy used for the M. capixaba ND2 case: search missing reference genes in long reads, combine gene-positive reads with a controlled set of mitogenome reads, and prepare a local assembly pool.\n")
    lines.append("## Read pools\n")
    for row in pool_rows:
        lines.append(
            f"- `{row['read_set']}`: mito_reads={row['mitogenome_read_ids_after_subsampling']} "
            f"(estimated_target={row.get('mitogenome_target_reads_requested', '.')}; "
            f"passing_filters={row.get('mitogenome_read_ids_passing_filters', '.')}) "
            f"+ missing_gene_reads={row['missing_gene_read_ids']} "
            f"=> pool={row['written_reads']} reads; approx_cov={row['approx_pool_coverage']}x; "
            f"filters: mito_aln_len>={row.get('mito_min_alignment_length', '.')}, "
            f"gene_aln_fraction>={row.get('gene_min_target_aligned_fraction', '.')}, "
            f"MAPQ>={row.get('min_mapq', '.')}; FASTQ={row['output_fastq']}; "
            f"flye_status={row.get('flye_status', '.')}; assembly={row.get('flye_assembly', '.')}"
        )
    (outdir / "missing_gene_recovery_report.md").write_text("\n".join(lines), encoding="utf-8")

    return outdir
