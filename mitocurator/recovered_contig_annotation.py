from __future__ import annotations

import csv
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

from .utils import ensure_dir, safe_get


def _write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _run_shell(cmd: str, log_path: Path, cwd: Optional[Path] = None) -> int:
    ensure_dir(log_path.parent)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(cmd + "\n\n")
        proc = subprocess.run(
            cmd,
            shell=True,
            cwd=str(cwd) if cwd else None,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return int(proc.returncode)


def _feature_name(feature: SeqFeature) -> str:
    q = feature.qualifiers
    for key in ["gene", "product", "locus_tag", "label"]:
        if key in q and q[key]:
            return str(q[key][0])
    return feature.type


def _feature_bounds(feature: SeqFeature) -> tuple[int, int]:
    loc = feature.location
    parts = getattr(loc, "parts", None)
    if parts:
        return min(int(p.start) for p in parts) + 1, max(int(p.end) for p in parts)
    return int(loc.start) + 1, int(loc.end)


def _safe_label(value: str) -> str:
    value = str(value or "sample").strip()
    out = []
    for char in value:
        if char.isalnum():
            out.append(char)
        elif char in {" ", "-", ".", "_"}:
            out.append("_")
    label = "".join(out).strip("_")
    while "__" in label:
        label = label.replace("__", "_")
    return label or "sample"


def _sample_label_from_config(config: dict) -> str:
    for path in [
        ["sample"],
        ["sample_name"],
        ["organism"],
        ["project", "sample"],
        ["project", "sample_name"],
        ["project", "organism"],
    ]:
        value = safe_get(config, path, None)
        if value:
            return _safe_label(str(value))
    return "sample"


def _normalize_selected_contig(in_fasta: Path, out_fasta: Path, sample: str) -> tuple[str, int]:
    if not in_fasta.exists():
        raise FileNotFoundError(f"Missing selected recovery contig FASTA: {in_fasta}")

    record = SeqIO.read(str(in_fasta), "fasta")
    old_id = record.id

    record = deepcopy(record)
    record.id = f"{sample}_recovered_contig"
    record.name = record.id
    record.description = f"{record.id} selected by MitoCurator missing-gene assembly assessment"

    ensure_dir(out_fasta.parent)
    SeqIO.write(record, str(out_fasta), "fasta")

    return old_id, len(record.seq)


def _mitofinder_command(config: dict, fasta: Path, job_name: str, workdir: Path) -> str:
    mitofinder = safe_get(config, ["mitofinder"], {}) or {}

    mode = str(mitofinder.get("mode", "python_interpreter"))
    script = str(mitofinder.get("script", "mitofinder"))
    python2 = str(mitofinder.get("python2", "python2"))
    reference_gb = str(mitofinder.get("reference_gb", ""))
    organism_code = str(mitofinder.get("organism_code", 5))
    threads = str(mitofinder.get("threads", safe_get(config, ["threads", "mitofinder"], 8)))

    if not reference_gb:
        raise ValueError("mitofinder.reference_gb is required for recovered contig annotation")

    base = script
    if mode == "python_interpreter":
        base = f"{python2} {script}"
    elif mode == "wrapper":
        base = str(mitofinder.get("wrapper", script))
    elif mode == "conda_env":
        conda = str(mitofinder.get("conda_executable", "conda"))
        env = str(mitofinder.get("conda_env", "mitofinder_py2"))
        base = f"{conda} run -n {env} {script}"

    # MitoFinder writes job-specific outputs inside the current working directory.
    # These are the standard options used by the pipeline: job name, assembly FASTA,
    # reference GenBank, mitochondrial genetic code and threads.
    return (
        f"{base} "
        f"-j {job_name} "
        f"-a {fasta} "
        f"-r {reference_gb} "
        f"-o {organism_code} "
        f"-p {threads}"
    )


def _find_best_genbank(workdir: Path) -> Optional[Path]:
    candidates = []
    for suffix in ("*.gb", "*.gbk", "*.genbank"):
        candidates.extend(workdir.rglob(suffix))

    parsed = []
    for path in candidates:
        try:
            rec = SeqIO.read(str(path), "genbank")
            cds_count = sum(1 for f in rec.features if f.type == "CDS")
            parsed.append((cds_count, len(rec.seq), path))
        except Exception:
            continue

    if not parsed:
        return None

    parsed.sort(key=lambda x: (x[0], x[1]), reverse=True)
    return parsed[0][2]


def _gene_inventory(record) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []

    for idx, feature in enumerate(record.features, start=1):
        if feature.type not in {"gene", "CDS", "rRNA", "tRNA", "misc_feature"}:
            continue

        start, end = _feature_bounds(feature)
        rows.append({
            "feature_index": idx,
            "type": feature.type,
            "gene": _feature_name(feature),
            "start": start,
            "end": end,
            "strand": int(feature.location.strand or 0),
            "length_nt": end - start + 1,
            "product": ";".join(feature.qualifiers.get("product", ["."])),
            "note": ";".join(feature.qualifiers.get("note", ["."])),
        })

    return rows


def _cds_validation(record, genetic_code: int) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene = _feature_name(feature)
        start, end = _feature_bounds(feature)
        strand = int(feature.location.strand or 0)
        nt = feature.extract(record.seq)
        aa = nt.translate(table=genetic_code, to_stop=False)

        aa_str = str(aa)
        internal = aa_str[:-1].count("*") if aa_str else 0
        internal_positions = [
            str(i + 1)
            for i, char in enumerate(aa_str[:-1])
            if char == "*"
        ]

        terminal_stop = "yes" if aa_str.endswith("*") else "no"
        length_multiple = "yes" if len(nt) % 3 == 0 else "no"

        if length_multiple == "yes" and internal == 0:
            status = "OK"
        elif internal > 0:
            status = "INTERNAL_STOP"
        else:
            status = "LENGTH_NOT_MULTIPLE_OF_3"

        rows.append({
            "gene": gene,
            "start": start,
            "end": end,
            "strand": strand,
            "cds_length_nt": len(nt),
            "length_multiple_of_3": length_multiple,
            "aa_length": len(aa_str),
            "terminal_stop": terminal_stop,
            "internal_stop_count": internal,
            "internal_stop_positions": ";".join(internal_positions) if internal_positions else ".",
            "status": status,
        })

    return rows


def _write_report(
    path: Path,
    selected_input: Path,
    renamed_fasta: Path,
    mitofinder_status: str,
    genbank: Optional[Path],
    inventory: List[Dict[str, Any]],
    cds_rows: List[Dict[str, Any]],
) -> None:
    cds = [r for r in inventory if r["type"] == "CDS"]
    rrna = [r for r in inventory if r["type"] == "rRNA"]
    trna = [r for r in inventory if r["type"] == "tRNA"]
    bad_cds = [r for r in cds_rows if r["status"] != "OK"]

    lines = []
    lines.append("# MitoCurator recovered contig annotation report\n")
    lines.append(f"- Selected contig input: `{selected_input}`")
    lines.append(f"- Renamed FASTA: `{renamed_fasta}`")
    lines.append(f"- MitoFinder status: {mitofinder_status}")
    lines.append(f"- GenBank selected: `{genbank}`" if genbank else "- GenBank selected: not found")
    lines.append("")
    lines.append("## Annotation summary\n")
    lines.append(f"- CDS: {len(cds)}")
    lines.append(f"- rRNAs: {len(rrna)}")
    lines.append(f"- tRNAs: {len(trna)}")
    lines.append(f"- CDS with remaining issues: {len(bad_cds)}")
    lines.append("")

    lines.append("## CDS validation\n")
    if cds_rows:
        for row in cds_rows:
            lines.append(
                f"- `{row['gene']}`: {row['status']}; "
                f"len={row['cds_length_nt']} nt; "
                f"internal_stops={row['internal_stop_count']}; "
                f"positions={row['internal_stop_positions']}"
            )
    else:
        lines.append("- No CDS validation rows generated.")

    path.write_text("\n".join(lines), encoding="utf-8")


def run_recovered_contig_annotation(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "11_recovered_contig_annotation")
    logs = ensure_dir(outdir / "logs")

    selected_contig = root / "10_missing_gene_assembly_assessment" / "selected_recovery_contig.fasta"
    sample = _sample_label_from_config(config)
    genetic_code = int(safe_get(config, ["mitofinder", "organism_code"], 5))

    renamed_fasta = outdir / f"{sample}_recovered_contig.fasta"
    old_id, contig_len = _normalize_selected_contig(selected_contig, renamed_fasta, sample)

    job_name = f"{sample}_recovered_contig"
    mf_workdir = ensure_dir(outdir / "mitofinder")
    cmd = _mitofinder_command(config, renamed_fasta, job_name, mf_workdir)

    (outdir / "run_mitofinder_recovered_contig.sh").write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"cd {mf_workdir}\n"
        f"{cmd}\n",
        encoding="utf-8",
    )
    (outdir / "run_mitofinder_recovered_contig.sh").chmod(0o755)

    run_mf = bool(safe_get(config, ["recovered_contig_annotation", "run_mitofinder"], True))
    returncode = None
    mitofinder_status = "not_run"

    if run_mf:
        returncode = _run_shell(cmd, logs / "mitofinder_recovered_contig.log", cwd=mf_workdir)
        mitofinder_status = "completed" if returncode == 0 else f"failed_returncode_{returncode}"

    genbank = _find_best_genbank(mf_workdir)

    inventory: List[Dict[str, Any]] = []
    cds_rows: List[Dict[str, Any]] = []

    if genbank is not None:
        record = SeqIO.read(str(genbank), "genbank")
        inventory = _gene_inventory(record)
        cds_rows = _cds_validation(record, genetic_code)

        # Keep a stable copy of the selected recovered annotation.
        SeqIO.write(record, str(outdir / "recovered_contig_annotated.gb"), "genbank")
        SeqIO.write(record, str(outdir / "recovered_contig_annotated.fasta"), "fasta")

    _write_tsv(
        outdir / "recovered_contig_annotation_inventory.tsv",
        inventory,
        [
            "feature_index",
            "type",
            "gene",
            "start",
            "end",
            "strand",
            "length_nt",
            "product",
            "note",
        ],
    )

    _write_tsv(
        outdir / "recovered_contig_cds_validation.tsv",
        cds_rows,
        [
            "gene",
            "start",
            "end",
            "strand",
            "cds_length_nt",
            "length_multiple_of_3",
            "aa_length",
            "terminal_stop",
            "internal_stop_count",
            "internal_stop_positions",
            "status",
        ],
    )

    _write_tsv(
        outdir / "recovered_contig_annotation_summary.tsv",
        [{
            "input_contig": str(selected_contig),
            "old_contig_id": old_id,
            "renamed_fasta": str(renamed_fasta),
            "contig_length": contig_len,
            "mitofinder_status": mitofinder_status,
            "mitofinder_returncode": returncode if returncode is not None else ".",
            "selected_genbank": str(genbank) if genbank else ".",
            "cds_count": sum(1 for r in inventory if r["type"] == "CDS"),
            "rrna_count": sum(1 for r in inventory if r["type"] == "rRNA"),
            "trna_count": sum(1 for r in inventory if r["type"] == "tRNA"),
            "cds_with_issues": sum(1 for r in cds_rows if r["status"] != "OK"),
        }],
        [
            "input_contig",
            "old_contig_id",
            "renamed_fasta",
            "contig_length",
            "mitofinder_status",
            "mitofinder_returncode",
            "selected_genbank",
            "cds_count",
            "rrna_count",
            "trna_count",
            "cds_with_issues",
        ],
    )

    _write_report(
        outdir / "recovered_contig_annotation_report.md",
        selected_contig,
        renamed_fasta,
        mitofinder_status,
        genbank,
        inventory,
        cds_rows,
    )

    return outdir
