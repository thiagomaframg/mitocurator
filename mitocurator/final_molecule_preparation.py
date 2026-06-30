from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from .io import infer_format
from .utils import ensure_dir, safe_get


_VALID_GB_DIVISIONS = {
    "BCT", "CON", "ENV", "EST", "GSS", "HTC", "HTG", "INV", "MAM", "PAT",
    "PHG", "PLN", "PRI", "ROD", "STS", "SYN", "TSA", "UNA", "VRL", "VRT",
}


def detect_at_rich_region(
    seq: str,
    window_size: int,
    min_at_pct: float,
    min_region_len: int,
) -> Optional[Dict[str, Any]]:
    """Find the most AT-rich contiguous region using a circular sliding window."""
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return None
    window_size = max(1, min(window_size, n))
    if min_region_len <= 0:
        min_region_len = window_size

    circular = seq + seq[: window_size - 1]
    mask = [False] * n

    at_count = sum(1 for b in circular[:window_size] if b in {"A", "T"})
    if (at_count / window_size) * 100 >= min_at_pct:
        mask[0] = True
    for i in range(1, n):
        if circular[i - 1] in {"A", "T"}:
            at_count -= 1
        if circular[i + window_size - 1] in {"A", "T"}:
            at_count += 1
        if (at_count / window_size) * 100 >= min_at_pct:
            mask[i] = True

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
        intervals.append((start, n - 1))
    if (len(intervals) > 1
            and intervals[0][0] == 0
            and intervals[-1][1] == n - 1):
        first = intervals.pop(0)
        last = intervals.pop(-1)
        intervals.insert(0, (last[0], first[1] + n))

    candidates: List[Dict[str, Any]] = []
    for s, e in intervals:
        length = e - s + 1
        if length < min_region_len:
            continue
        region_seq = (seq + seq)[s : e + 1]
        at_pct = sum(1 for b in region_seq if b in {"A", "T"}) / len(region_seq) * 100
        candidates.append({
            "start": (s % n) + 1,
            "end": (e % n) + 1,
            "length": length,
            "at_percent": round(at_pct, 3),
            "wraps_around_origin": "yes" if e >= n else "no",
            "_score": (at_pct, length),
        })

    if not candidates:
        return None
    candidates.sort(key=lambda x: x["_score"], reverse=True)
    best = candidates[0]
    del best["_score"]
    return best


def _feature_inventory(record) -> List[Dict[str, Any]]:
    rows = []
    for feat in record.features:
        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        strand = int(feat.location.strand or 0)
        gene = ";".join(feat.qualifiers.get("gene", [])) or "."
        product = ";".join(feat.qualifiers.get("product", [])) or "."
        notes = ";".join(feat.qualifiers.get("note", [])) or "."
        rows.append({
            "feature_type": feat.type,
            "gene": gene,
            "product": product,
            "start": start,
            "end": end,
            "strand": strand,
            "length_nt": end - start + 1,
            "qualifiers_summary": notes,
        })
    return rows


def _write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_final_molecule_preparation(
    config: dict,
    input_gb: Path,
    outdir: Path,
    organism: Optional[str] = None,
    topology: Optional[str] = None,
    expected_length: Optional[int] = None,
    at_window_size: int = 500,
    at_min_pct: float = 85.0,
    at_min_len: int = 500,
) -> Path:
    ensure_dir(outdir)

    # genetic_code: required — silent fallback would produce wrong translations for non-invertebrates
    tt = safe_get(config, ["project", "genetic_code"], None)
    if tt is None:
        raise KeyError(
            "project.genetic_code is required in config.yaml "
            "(e.g. genetic_code: 5 for invertebrate mitochondria, 2 for vertebrate)"
        )
    tt = int(tt)

    # gb_division: required — INV/VRT/PLN differ and affect the GenBank LOCUS line
    gb_division = safe_get(config, ["project", "gb_division"], None)
    if gb_division is None:
        raise KeyError(
            "project.gb_division is required in config.yaml "
            f"(e.g. gb_division: INV for invertebrates, VRT for vertebrates). "
            f"Known values: {', '.join(sorted(_VALID_GB_DIVISIONS))}"
        )
    gb_division = str(gb_division).upper()
    if gb_division not in _VALID_GB_DIVISIONS:
        raise ValueError(
            f"project.gb_division {gb_division!r} is not a recognised GenBank division. "
            f"Known values: {', '.join(sorted(_VALID_GB_DIVISIONS))}"
        )

    if not input_gb.exists():
        raise FileNotFoundError(f"Input GenBank not found: {input_gb}")
    fmt = infer_format(input_gb)
    record = SeqIO.read(str(input_gb), fmt)

    if topology is None:
        topology = str(safe_get(config, ["project", "topology"], "circular"))

    if not organism:
        for keys in (["project", "organism"], ["organism"], ["sample_name"]):
            v = safe_get(config, keys, None)
            if v:
                organism = str(v)
                break
    if not organism:
        organism = "Unknown organism"

    # 2. Set GenBank LOCUS metadata
    record.annotations["topology"] = topology
    record.annotations["molecule_type"] = "DNA"
    record.annotations["data_file_division"] = gb_division
    record.annotations.setdefault("keywords", ["."])

    # 3. Recreate source feature spanning the full sequence
    source_feat = SeqFeature(
        location=FeatureLocation(0, len(record.seq), strand=1),
        type="source",
        qualifiers={
            "organism": [organism],
            "organelle": ["mitochondrion"],
            "mol_type": ["genomic DNA"],
        },
    )

    # 4. Inject transl_table into all CDS features
    non_source = [f for f in record.features if f.type != "source"]
    for feat in non_source:
        if feat.type == "CDS":
            feat.qualifiers["transl_table"] = [str(tt)]

    # 5. Detect A+T-rich region and annotate as misc_feature
    at_region = detect_at_rich_region(str(record.seq), at_window_size, at_min_pct, at_min_len)
    if at_region:
        s0 = at_region["start"] - 1
        e0 = at_region["end"]
        if at_region["wraps_around_origin"] == "yes":
            loc = (FeatureLocation(s0, len(record.seq), strand=1)
                   + FeatureLocation(0, e0, strand=1))
        else:
            loc = FeatureLocation(s0, e0, strand=1)
        at_feat = SeqFeature(
            location=loc,
            type="misc_feature",
            qualifiers={"note": [
                "A+T-rich control region",
                "putative mitochondrial control region",
            ]},
        )
        non_source.append(at_feat)

    record.features = [source_feat] + non_source

    # 6. Write output files
    fasta_out  = outdir / "final_molecule.fasta"
    gb_out     = outdir / "final_molecule.gb"
    inv_out    = outdir / "final_molecule_feature_inventory.tsv"
    at_out     = outdir / "at_rich_region.tsv"
    report_out = outdir / "final_molecule_report.md"

    SeqIO.write(record, str(fasta_out), "fasta")
    SeqIO.write(record, str(gb_out), "genbank")

    inventory = _feature_inventory(record)
    _write_tsv(
        inv_out, inventory,
        ["feature_type", "gene", "product", "start", "end", "strand", "length_nt", "qualifiers_summary"],
    )

    if at_region:
        at_rows: List[Dict[str, Any]] = [at_region]
        at_fields = list(at_region.keys())
    else:
        at_rows = [{"start": ".", "end": ".", "length": 0, "at_percent": ".",
                    "wraps_around_origin": "no", "status": "not detected"}]
        at_fields = list(at_rows[0].keys())
    _write_tsv(at_out, at_rows, at_fields)

    # 7. Markdown report
    cds_count  = sum(1 for f in record.features if f.type == "CDS")
    rrna_count = sum(1 for f in record.features if f.type == "rRNA")
    trna_count = sum(1 for f in record.features if f.type == "tRNA")

    lines = [
        "# MitoCurator final molecule preparation report\n",
        f"- Input: `{input_gb}` ({fmt})",
        f"- Organism: {organism}",
        f"- Final molecule length: {len(record.seq)} bp",
        f"- Topology: {topology}",
        f"- GenBank division: {gb_division}",
        f"- Genetic code (transl_table): {tt}",
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
            f"- A+T-rich/control region: {at_region['start']}..{at_region['end']} "
            f"(length={at_region['length']} bp, AT={at_region['at_percent']}%)",
            f"- Wraps around origin: {at_region['wraps_around_origin']}",
        ])
    else:
        lines.append("- A+T-rich/control region: not detected with configured thresholds")

    lines.extend([
        "",
        "## Interpretation guardrails\n",
        "- This step prepares a candidate final mitochondrial molecule with explicit circular mitochondrial genomic DNA annotation.",
        "- Biological circularity and absence of CNVs are not proven here; coverage and CNV validation are handled separately.",
    ])

    report_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return gb_out
