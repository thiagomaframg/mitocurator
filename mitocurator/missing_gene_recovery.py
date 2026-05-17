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
    target_query_bases: int,
    min_reads: int,
    rng: random.Random,
) -> Set[str]:
    """Select mitogenome reads for a controlled local assembly pool.

    The M. capixaba manual workflow used a practical pool size, not a
    requirement that every read spans a large fraction of the mitogenome.
    Therefore, this selector accumulates full read bases and also enforces
    a minimum number of mitogenome-like reads.
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

    # Prefer longer/better alignments, but shuffle within the candidate set to
    # avoid always taking one local overrepresented molecule first.
    rng.shuffle(rows)
    rows.sort(
        key=lambda r: (
            int(r["alignment_length"]),
            float(r["identity"]),
            int(r["mapq"]),
        ),
        reverse=True,
    )

    selected: Set[str] = set()
    query_bases = 0

    for row in rows:
        rid = str(row["read_id"])
        if rid in selected:
            continue

        selected.add(rid)
        query_bases += int(row["query_length"])

        if len(selected) >= min_reads and query_bases >= target_query_bases:
            break

    return selected


def _select_read_ids_from_paf(
    paf_rows: List[Dict[str, Any]],
    min_mapq: int,
    min_identity: float,
    min_target_coverage: float,
) -> Set[str]:
    keep = set()
    for row in paf_rows:
        if int(row["mapq"]) < min_mapq:
            continue
        if float(row["identity"]) < min_identity:
            continue
        if float(row["target_coverage"]) < min_target_coverage:
            continue
        keep.add(str(row["read_id"]))
    return keep


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
    min_mapq = int(safe_get(config, ["missing_gene_recovery", "min_mapq"], 20))
    gene_min_identity = float(safe_get(config, ["missing_gene_recovery", "gene_min_identity"], 70))
    gene_min_target_cov = float(safe_get(config, ["missing_gene_recovery", "gene_min_target_coverage"], 50))
    mito_min_identity = float(safe_get(config, ["missing_gene_recovery", "mito_min_identity"], 70))
    mito_min_alignment_length = int(safe_get(config, ["missing_gene_recovery", "mito_min_alignment_length"], 500))
    mito_min_reads = int(safe_get(config, ["missing_gene_recovery", "mito_min_reads"], 300))
    seed = int(safe_get(config, ["missing_gene_recovery", "random_seed"], 42))
    rng = random.Random(seed)

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

        gene_ids = _select_read_ids_from_paf(
            gene_paf_rows,
            min_mapq=min_mapq,
            min_identity=gene_min_identity,
            min_target_coverage=gene_min_target_cov,
        )
        target_query_bases = int(target_cov * mito_len)
        mito_ids_all = _select_mitogenome_read_ids_for_target_coverage(
            mito_paf_rows,
            min_mapq=min_mapq,
            min_identity=mito_min_identity,
            min_alignment_length=mito_min_alignment_length,
            target_query_bases=target_query_bases,
            min_reads=mito_min_reads,
            rng=rng,
        )

        n_reads, n_bases = _fastq_stats(rs["reads"])
        mean_read_len = n_bases / n_reads if n_reads else 0

        mito_ids = sorted(mito_ids_all)
        final_ids = set(mito_ids) | gene_ids

        out_fastq = selected_dir / f"{name}.missing_gene_recovery_pool.fastq.gz"
        written_reads, written_bases = _write_reads_by_id(rs["reads"], out_fastq, final_ids)

        flye_script = selected_dir / f"{name}.missing_gene_recovery.flye.sh"
        flye_script.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n\n"
            f"flye --pacbio-hifi {out_fastq} --genome-size 20k "
            f"--threads {threads} --out-dir {selected_dir / (name + '.flye_missing_gene_recovery')}\n",
            encoding="utf-8",
        )
        flye_script.chmod(0o755)

        pool_rows.append({
            "read_set": name,
            "read_type": rtype,
            "input_reads": n_reads,
            "input_bases": n_bases,
            "mitogenome_length": mito_len,
            "target_mitogenome_coverage": target_cov,
            "mean_read_length": round(mean_read_len, 3),
            "mitogenome_read_ids_before_subsampling": ".",
            "mitogenome_read_ids_after_subsampling": len(mito_ids),
            "missing_gene_read_ids": len(gene_ids),
            "final_pool_read_ids": len(final_ids),
            "written_reads": written_reads,
            "written_bases": written_bases,
            "approx_pool_coverage": round(written_bases / mito_len, 3) if mito_len else 0,
            "output_fastq": str(out_fastq),
            "flye_script": str(flye_script),
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
            "mitogenome_read_ids_before_subsampling",
            "mitogenome_read_ids_after_subsampling",
            "missing_gene_read_ids",
            "final_pool_read_ids",
            "written_reads",
            "written_bases",
            "approx_pool_coverage",
            "output_fastq",
            "flye_script",
        ],
    )

    lines = []
    lines.append("# MitoCurator missing-gene recovery report\n")
    lines.append("This step follows the direct recovery strategy used for the M. capixaba ND2 case: search missing reference genes in long reads, combine gene-positive reads with a controlled set of mitogenome reads, and prepare a local assembly pool.\n")
    lines.append("## Read pools\n")
    for row in pool_rows:
        lines.append(
            f"- `{row['read_set']}`: mito_reads={row['mitogenome_read_ids_after_subsampling']} "
            f"+ missing_gene_reads={row['missing_gene_read_ids']} "
            f"=> pool={row['written_reads']} reads; approx_cov={row['approx_pool_coverage']}x; "
            f"FASTQ={row['output_fastq']}"
        )
    (outdir / "missing_gene_recovery_report.md").write_text("\n".join(lines), encoding="utf-8")

    return outdir
