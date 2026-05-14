from __future__ import annotations

import csv
import gzip
import subprocess
from pathlib import Path
from typing import Any, Dict, List

from Bio import SeqIO

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


def _safe_read_tsv(path: Path) -> List[Dict[str, str]]:
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


def _feature_name(feature) -> str:
    q = feature.qualifiers
    for key in ["gene", "product", "locus_tag", "label"]:
        if key in q and q[key]:
            return str(q[key][0])
    return feature.type


def _feature_positions(feature) -> set[int]:
    positions: set[int] = set()
    parts = getattr(feature.location, "parts", [feature.location])
    for part in parts:
        start = int(part.start) + 1
        end = int(part.end)
        positions.update(range(start, end + 1))
    return positions


def _read_features(root: Path) -> List[Dict[str, Any]]:
    gb = root / "05_refinement" / "refined.gb"
    if not gb.exists():
        return []

    record = SeqIO.read(str(gb), "genbank")
    features: List[Dict[str, Any]] = []

    for feature in record.features:
        if feature.type not in {"CDS", "rRNA", "tRNA", "misc_feature"}:
            continue
        positions = _feature_positions(feature)
        if not positions:
            continue

        features.append({
            "gene": _feature_name(feature),
            "type": feature.type,
            "start": min(positions),
            "end": max(positions),
            "positions": positions,
        })

    return features


def _locate_feature(pos: int, features: List[Dict[str, Any]]) -> Dict[str, str]:
    hits = [f for f in features if pos in f["positions"]]
    if not hits:
        return {
            "gene": ".",
            "feature_type": "intergenic",
            "feature_start": ".",
            "feature_end": ".",
        }

    # Prefer CDS/rRNA/tRNA over misc_feature when overlapping.
    priority = {"CDS": 0, "rRNA": 1, "tRNA": 2, "misc_feature": 3}
    hits.sort(key=lambda f: priority.get(f["type"], 99))
    hit = hits[0]

    return {
        "gene": hit["gene"],
        "feature_type": hit["type"],
        "feature_start": str(hit["start"]),
        "feature_end": str(hit["end"]),
    }


def _variant_type(ref: str, alt: str) -> str:
    alts = alt.split(",")
    if any(len(a) != len(ref) for a in alts):
        return "INDEL"
    if len(ref) == 1 and all(len(a) == 1 for a in alts):
        return "SNP"
    return "MNP"


def _parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "true"
    return out


def _parse_sample(format_field: str, sample_field: str) -> Dict[str, str]:
    keys = format_field.split(":")
    values = sample_field.split(":")
    return {k: v for k, v in zip(keys, values)}


def _alt_depth_from_ad(ad: str) -> str:
    if not ad or ad == ".":
        return "."
    parts = ad.split(",")
    if len(parts) < 2:
        return "."
    try:
        return str(sum(int(x) for x in parts[1:] if x != "."))
    except ValueError:
        return "."


def _alt_freq(dp: str, ad: str) -> str:
    alt_dp = _alt_depth_from_ad(ad)
    if alt_dp == "." or not dp or dp == ".":
        return "."
    try:
        dp_i = int(dp)
        alt_i = int(alt_dp)
    except ValueError:
        return "."
    if dp_i == 0:
        return "."
    return f"{alt_i / dp_i:.4f}"


def _read_vcf(vcf_gz: Path, readset: str, features: List[Dict[str, Any]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []

    with gzip.open(vcf_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _id, ref, alt, qual, filt, info = parts[:8]
            fmt = parts[8] if len(parts) > 8 else ""
            sample = parts[9] if len(parts) > 9 else ""

            info_dict = _parse_info(info)
            sample_dict = _parse_sample(fmt, sample) if fmt and sample else {}

            dp = sample_dict.get("DP") or info_dict.get("DP", ".")
            ad = sample_dict.get("AD", ".")
            alt_dp = _alt_depth_from_ad(ad)
            af = _alt_freq(dp, ad)

            feature = _locate_feature(int(pos), features)
            vtype = _variant_type(ref, alt)

            rows.append({
                "readset": readset,
                "chrom": chrom,
                "pos": pos,
                "gene": feature["gene"],
                "feature_type": feature["feature_type"],
                "feature_start": feature["feature_start"],
                "feature_end": feature["feature_end"],
                "ref": ref,
                "alt": alt,
                "variant_type": vtype,
                "qual": qual,
                "filter": filt,
                "depth": dp,
                "allelic_depth": ad,
                "alt_depth": alt_dp,
                "alt_frequency": af,
            })

    return rows


def _summarize_gene_variants(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    grouped: Dict[tuple[str, str, str], List[Dict[str, str]]] = {}

    for row in rows:
        key = (row["readset"], row["gene"], row["feature_type"])
        grouped.setdefault(key, []).append(row)

    out: List[Dict[str, str]] = []
    for (readset, gene, feature_type), items in sorted(grouped.items()):
        n_total = len(items)
        n_snp = sum(1 for r in items if r["variant_type"] == "SNP")
        n_indel = sum(1 for r in items if r["variant_type"] == "INDEL")
        max_af = "."
        afs = []
        for r in items:
            try:
                afs.append(float(r["alt_frequency"]))
            except ValueError:
                pass
        if afs:
            max_af = f"{max(afs):.4f}"

        high_af_positions = [
            r["pos"] for r in items
            if r["alt_frequency"] != "." and float(r["alt_frequency"]) >= 0.8
        ]

        out.append({
            "readset": readset,
            "gene": gene,
            "feature_type": feature_type,
            "n_variants": str(n_total),
            "n_snps": str(n_snp),
            "n_indels": str(n_indel),
            "max_alt_frequency": max_af,
            "high_alt_frequency_positions": ";".join(high_af_positions) if high_af_positions else ".",
        })

    return out


def run_variant_evidence(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "09_variant_evidence")
    logs = ensure_dir(outdir / "logs")

    read_mapping_dir = root / "08_read_mapping"
    ref_fasta = read_mapping_dir / "mitogenome_reference.fasta"
    readsets_tsv = read_mapping_dir / "readsets.tsv"

    if not ref_fasta.exists():
        raise FileNotFoundError(f"Missing reference FASTA from read mapping stage: {ref_fasta}")
    if not readsets_tsv.exists():
        raise FileNotFoundError(f"Missing readsets.tsv from read mapping stage: {readsets_tsv}")

    samtools = str(safe_get(config, ["tools", "samtools"], "samtools"))
    bcftools = str(safe_get(config, ["tools", "bcftools"], "bcftools"))

    threads = int(safe_get(config, ["variant_evidence", "threads"], safe_get(config, ["threads", "bcftools"], 4)))
    min_depth = int(safe_get(config, ["variant_evidence", "min_depth"], safe_get(config, ["targeted_consensus", "min_depth"], 5)))
    min_qual = float(safe_get(config, ["variant_evidence", "min_qual"], 0))

    _run_shell(f"{samtools} faidx {ref_fasta}", logs / "reference.faidx.log")

    readsets = _safe_read_tsv(readsets_tsv)
    features = _read_features(root)

    all_variant_rows: List[Dict[str, str]] = []

    for rs in readsets:
        name = rs["name"]
        bam = read_mapping_dir / f"{name}.sorted.bam"

        if not bam.exists():
            raise FileNotFoundError(f"Missing BAM for readset '{name}': {bam}")

        raw_vcf = outdir / f"{name}.raw.vcf.gz"
        filtered_vcf = outdir / f"{name}.filtered.vcf.gz"

        cmd_raw = (
            f"{bcftools} mpileup --threads {threads} "
            f"-a FORMAT/DP,FORMAT/AD -Ou -f {ref_fasta} {bam} "
            f"| {bcftools} call --threads {threads} -mv -Oz -o {raw_vcf}"
        )
        _run_shell(cmd_raw, logs / f"{name}.mpileup_call.log")

        _run_shell(f"{bcftools} index -t {raw_vcf}", logs / f"{name}.raw_index.log")

        filter_expr = f"INFO/DP>={min_depth} && QUAL>={min_qual}"
        cmd_filter = f"{bcftools} view -i '{filter_expr}' -Oz -o {filtered_vcf} {raw_vcf}"
        _run_shell(cmd_filter, logs / f"{name}.filter.log")

        _run_shell(f"{bcftools} index -t {filtered_vcf}", logs / f"{name}.filtered_index.log")

        rows = _read_vcf(filtered_vcf, name, features)

        _write_tsv(
            outdir / f"{name}.variant_summary.tsv",
            rows,
            [
                "readset",
                "chrom",
                "pos",
                "gene",
                "feature_type",
                "feature_start",
                "feature_end",
                "ref",
                "alt",
                "variant_type",
                "qual",
                "filter",
                "depth",
                "allelic_depth",
                "alt_depth",
                "alt_frequency",
            ],
        )

        all_variant_rows.extend(rows)

    _write_tsv(
        outdir / "variant_summary.tsv",
        all_variant_rows,
        [
            "readset",
            "chrom",
            "pos",
            "gene",
            "feature_type",
            "feature_start",
            "feature_end",
            "ref",
            "alt",
            "variant_type",
            "qual",
            "filter",
            "depth",
            "allelic_depth",
            "alt_depth",
            "alt_frequency",
        ],
    )

    gene_rows = _summarize_gene_variants(all_variant_rows)
    _write_tsv(
        outdir / "gene_variant_evidence.tsv",
        gene_rows,
        [
            "readset",
            "gene",
            "feature_type",
            "n_variants",
            "n_snps",
            "n_indels",
            "max_alt_frequency",
            "high_alt_frequency_positions",
        ],
    )

    return outdir
