from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from Bio import SeqIO

from .utils import ensure_dir, safe_get


def _write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


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


def _fasta_lengths(path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    if not path.exists():
        return rows

    for rec in SeqIO.parse(str(path), "fasta"):
        rows.append({
            "contig": rec.id,
            "length": len(rec.seq),
        })

    rows.sort(key=lambda r: int(r["length"]), reverse=True)
    return rows


def _parse_paf(path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []

    if not path.exists():
        return rows

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            # With command: minimap2 reference_genes.fasta assembly.fasta
            # PAF query = contig; target = missing reference gene.
            contig = parts[0]
            contig_len = int(parts[1])
            contig_start = int(parts[2])
            contig_end = int(parts[3])
            strand = parts[4]
            gene = parts[5]
            gene_len = int(parts[6])
            gene_start = int(parts[7])
            gene_end = int(parts[8])
            matches = int(parts[9])
            aln_len = int(parts[10])
            mapq = int(parts[11])

            gene_cov_pct = 100 * (gene_end - gene_start) / gene_len if gene_len else 0.0
            identity_pct = 100 * matches / aln_len if aln_len else 0.0

            rows.append({
                "contig": contig,
                "contig_len": contig_len,
                "contig_start": contig_start,
                "contig_end": contig_end,
                "strand": strand,
                "gene": gene,
                "gene_len": gene_len,
                "gene_start": gene_start,
                "gene_end": gene_end,
                "matches": matches,
                "aln_len": aln_len,
                "mapq": mapq,
                "gene_cov_pct": round(gene_cov_pct, 3),
                "identity_pct": round(identity_pct, 3),
            })

    return rows


def _decision(row: Dict[str, Any], config: dict) -> str:
    min_gene_cov = float(safe_get(config, ["missing_gene_assembly_assessment", "min_gene_coverage_pct"], 90.0))
    min_identity = float(safe_get(config, ["missing_gene_assembly_assessment", "min_identity_pct"], 65.0))
    min_mapq = int(safe_get(config, ["missing_gene_assembly_assessment", "min_mapq"], 30))

    min_contig_len = int(safe_get(config, ["missing_gene_assembly_assessment", "min_contig_length"], 15000))
    max_contig_len = int(safe_get(config, ["missing_gene_assembly_assessment", "max_contig_length"], 25000))

    if float(row["gene_cov_pct"]) < min_gene_cov:
        return "LOW_GENE_COVERAGE"
    if float(row["identity_pct"]) < min_identity:
        return "LOW_IDENTITY"
    if int(row["mapq"]) < min_mapq:
        return "LOW_MAPQ"

    contig_len = int(row["contig_len"])
    if min_contig_len <= contig_len <= max_contig_len:
        return "RECOVERY_CONTIG_SUPPORTED"

    return "GENE_SUPPORTED_BUT_CONTIG_SIZE_OUTSIDE_EXPECTED_RANGE"


def _score(row: Dict[str, Any], decision: str) -> float:
    score = 0.0

    if decision == "RECOVERY_CONTIG_SUPPORTED":
        score += 100000
    elif decision == "GENE_SUPPORTED_BUT_CONTIG_SIZE_OUTSIDE_EXPECTED_RANGE":
        score += 50000

    score += float(row["gene_cov_pct"]) * 100
    score += float(row["identity_pct"]) * 10
    score += int(row["mapq"])

    contig_len = int(row["contig_len"])
    # Prefer contigs close to 20 kb, but do not make this dominant over gene recovery.
    score -= abs(contig_len - 20000) / 1000

    return round(score, 3)


def _find_latest_assembly(root: Path, config: dict) -> Optional[Path]:
    explicit = safe_get(config, ["missing_gene_assembly_assessment", "assembly_fasta"], None)
    if explicit:
        p = Path(str(explicit))
        if p.exists():
            return p

    selected_reads = root / "09_missing_gene_recovery" / "selected_reads"
    if not selected_reads.exists():
        return None

    assemblies = sorted(
        selected_reads.glob("*.flye_missing_gene_recovery*/assembly.fasta"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )

    if assemblies:
        return assemblies[0]

    return None


def _extract_selected_contig(assembly: Path, contig_name: str, out_fasta: Path) -> bool:
    ensure_dir(out_fasta.parent)

    for rec in SeqIO.parse(str(assembly), "fasta"):
        if rec.id == contig_name:
            SeqIO.write(rec, str(out_fasta), "fasta")
            return True

    return False


def _write_report(path: Path, assembly: Optional[Path], rows: List[Dict[str, Any]], selected: Optional[Dict[str, Any]]) -> None:
    lines = []
    lines.append("# MitoCurator missing-gene assembly assessment\n")

    if assembly is None:
        lines.append("No Flye assembly was found for assessment.\n")
        path.write_text("\n".join(lines), encoding="utf-8")
        return

    lines.append(f"- Assembly evaluated: `{assembly}`")
    lines.append(f"- Candidate alignments: {len(rows)}\n")

    if selected:
        lines.append("## Selected recovery contig\n")
        lines.append(f"- Contig: `{selected['contig']}`")
        lines.append(f"- Length: {selected['contig_len']} bp")
        lines.append(f"- Gene: {selected['gene']}")
        lines.append(f"- Gene reference length: {selected['gene_len']} nt")
        lines.append(f"- Gene aligned length: {selected['aln_len']} nt")
        lines.append(f"- Gene coverage: {selected['gene_cov_pct']}%")
        lines.append(f"- Identity: {selected['identity_pct']}%")
        lines.append(f"- MAPQ: {selected['mapq']}")
        lines.append(f"- Strand: {selected['strand']}")
        lines.append(f"- Decision: {selected['decision']}")
        lines.append("")
    else:
        lines.append("## Selected recovery contig\n")
        lines.append("- No recovery contig passed the configured criteria.\n")

    lines.append("## Ranked alignments\n")
    for row in rows:
        lines.append(
            f"- `{row['contig']}` vs `{row['gene']}`: "
            f"gene_cov={row['gene_cov_pct']}%; "
            f"identity={row['identity_pct']}%; "
            f"mapq={row['mapq']}; "
            f"contig_len={row['contig_len']}; "
            f"decision={row['decision']}"
        )

    path.write_text("\n".join(lines), encoding="utf-8")


def run_missing_gene_assembly_assessment(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "10_missing_gene_assembly_assessment")
    logs = ensure_dir(outdir / "logs")

    ref_genes = root / "09_missing_gene_recovery" / "missing_reference_genes.fasta"
    if not ref_genes.exists():
        raise FileNotFoundError(f"Missing reference gene FASTA: {ref_genes}")

    assembly = _find_latest_assembly(root, config)

    contig_rows: List[Dict[str, Any]] = []
    aln_rows: List[Dict[str, Any]] = []
    selected: Optional[Dict[str, Any]] = None

    if assembly is not None and assembly.exists():
        contig_rows = _fasta_lengths(assembly)

        paf = outdir / "missing_gene_vs_recovery_assembly.paf"
        minimap2 = str(safe_get(config, ["tools", "minimap2"], "minimap2"))

        _run_shell(
            f"{minimap2} -x asm20 {ref_genes} {assembly} > {paf}",
            logs / "missing_gene_vs_recovery_assembly.log",
        )

        aln_rows = _parse_paf(paf)

        for row in aln_rows:
            decision = _decision(row, config)
            row["decision"] = decision
            row["score"] = _score(row, decision)

        aln_rows.sort(
            key=lambda r: (
                float(r["score"]),
                float(r["gene_cov_pct"]),
                float(r["identity_pct"]),
                -abs(int(r["contig_len"]) - 20000),
            ),
            reverse=True,
        )

        supported = [
            r for r in aln_rows
            if r["decision"] in {
                "RECOVERY_CONTIG_SUPPORTED",
                "GENE_SUPPORTED_BUT_CONTIG_SIZE_OUTSIDE_EXPECTED_RANGE",
            }
        ]

        if supported:
            selected = supported[0]
            _extract_selected_contig(
                assembly,
                str(selected["contig"]),
                outdir / "selected_recovery_contig.fasta",
            )

    _write_tsv(
        outdir / "flye_contig_stats.tsv",
        contig_rows,
        ["contig", "length"],
    )

    _write_tsv(
        outdir / "missing_gene_vs_recovery_assembly.tsv",
        aln_rows,
        [
            "contig",
            "contig_len",
            "contig_start",
            "contig_end",
            "strand",
            "gene",
            "gene_len",
            "gene_start",
            "gene_end",
            "matches",
            "aln_len",
            "mapq",
            "gene_cov_pct",
            "identity_pct",
            "decision",
            "score",
        ],
    )

    selected_rows = [selected] if selected else []
    _write_tsv(
        outdir / "selected_recovery_contig.tsv",
        selected_rows,
        [
            "contig",
            "contig_len",
            "contig_start",
            "contig_end",
            "strand",
            "gene",
            "gene_len",
            "gene_start",
            "gene_end",
            "matches",
            "aln_len",
            "mapq",
            "gene_cov_pct",
            "identity_pct",
            "decision",
            "score",
        ],
    )

    _write_report(
        outdir / "missing_gene_assembly_assessment_report.md",
        assembly,
        aln_rows,
        selected,
    )

    return outdir
