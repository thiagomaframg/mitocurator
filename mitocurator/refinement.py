from __future__ import annotations
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from .io import read_record, write_record
from .utils import ensure_dir, safe_get

EXPECTED_GENE_SET = {
    "CDS": ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"],
    "rRNA": ["rrnL", "rrnS"],
    "tRNA": [
        "tRNA-Ala", "tRNA-Arg", "tRNA-Asn", "tRNA-Asp", "tRNA-Cys", "tRNA-Gln", "tRNA-Glu", "tRNA-Gly", "tRNA-His",
        "tRNA-Ile", "tRNA-Leu", "tRNA-Leu2", "tRNA-Lys", "tRNA-Met", "tRNA-Phe", "tRNA-Pro", "tRNA-Ser", "tRNA-Ser2",
        "tRNA-Thr", "tRNA-Trp", "tRNA-Tyr", "tRNA-Val",
    ],
}


def _normalize_token(value: str) -> str:
    x = (value or "").strip()
    u = x.upper().replace("_", "").replace(" ", "")
    aliases = {
        "NAD1": "ND1", "NAD2": "ND2", "NAD3": "ND3", "NAD4": "ND4", "NAD4L": "ND4L", "NAD5": "ND5", "NAD6": "ND6",
        "CYTB": "CYTB", "COB": "CYTB", "COX1": "COX1", "COI": "COX1", "COXI": "COX1", "COX2": "COX2", "COII": "COX2", "COX3": "COX3", "COIII": "COX3",
        "RRNL": "rrnL", "RRNS": "rrnS",
        "TRNW": "tRNA-Trp",
        "TRNL1": "tRNA-Leu", "TRNL2": "tRNA-Leu2", "TRNS1": "tRNA-Ser", "TRNS2": "tRNA-Ser2",
    }
    if u in aliases:
        return aliases[u]
    if u.startswith("TRNA"):
        y = x.replace("_", "-").replace(" ", "").replace("--", "-")
        y = y.replace("tRNA", "tRNA-") if y.startswith("tRNA") and not y.startswith("tRNA-") else y
        y = y.replace("TRNA", "tRNA-") if y.startswith("TRNA") else y
        y = y.replace("-1", "1").replace("-2", "2")
        y = y.replace("Leu1", "Leu").replace("Leu-1", "Leu").replace("Leu2", "Leu2")
        y = y.replace("Ser1", "Ser").replace("Ser-1", "Ser").replace("Ser2", "Ser2")
        return y
    return x.upper()


def summarize_expected_gene_set(record, out_tsv: Path):
    counts = Counter()
    for feat in record.features:
        ftype = feat.type
        vals = []
        if ftype == "CDS":
            vals.extend(feat.qualifiers.get("gene", []))
        elif ftype in ("rRNA", "tRNA"):
            for k in ("gene", "product", "note"):
                vals.extend(feat.qualifiers.get(k, []))
        else:
            continue
        for v in vals:
            n = _normalize_token(v)
            counts[n] += 1

    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("gene\ttype\tstatus\tcopies\tcomment\n")
        for gtype, genes in EXPECTED_GENE_SET.items():
            for gene in genes:
                c = counts.get(gene, 0)
                status = "PRESENT" if c == 1 else "DUPLICATED" if c > 1 else "MISSING"
                comment = "." if c else "not detected in current annotation"
                out.write(f"{gene}\t{gtype}\t{status}\t{c}\t{comment}\n")


def _feature_segments(feat) -> List[Tuple[int, int]]:
    if isinstance(feat.location, CompoundLocation):
        return [(int(p.start), int(p.end)) for p in feat.location.parts]
    return [(int(feat.location.start), int(feat.location.end))]


def add_at_rich_region(record, min_len=500, min_at=75.0):
    n = len(record.seq)
    segments = []
    labels = []
    for feat in record.features:
        if feat.type == "source":
            continue
        segments.extend(_feature_segments(feat))
        labels.append(feat.qualifiers.get("gene", [feat.type])[0])
    if not segments:
        return None

    merged = []
    for s, e in sorted(segments):
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)

    gaps = []
    for i in range(len(merged) - 1):
        a, b = merged[i], merged[i + 1]
        if b[0] > a[1]:
            gaps.append((a[1], b[0], False))
    if merged[0][0] > merged[-1][1]:
        pass
    wrap_len = (merged[0][0] + n) - merged[-1][1]
    if wrap_len > 0:
        gaps.append((merged[-1][1], merged[0][0], True))
    if not gaps:
        return None

    best = max(gaps, key=lambda g: (g[1] - g[0]) if not g[2] else ((n - g[0]) + g[1]))
    s, e, wrap = best
    seq = (record.seq[s:e] if not wrap else (record.seq[s:] + record.seq[:e]))
    length = len(seq)
    at = 100.0 * (str(seq).upper().count("A") + str(seq).upper().count("T")) / length if length else 0.0
    if length < min_len or at < min_at:
        return None

    if wrap:
        loc = CompoundLocation([FeatureLocation(s, n, strand=1), FeatureLocation(0, e, strand=1)])
        start_1b, end_1b = s + 1, e
    else:
        loc = FeatureLocation(s, e, strand=1)
        start_1b, end_1b = s + 1, e

    feat = SeqFeature(
        location=loc,
        type="misc_feature",
        qualifiers={
            "note": ["putative AT-rich control region"],
            "product": ["AT-rich region"],
            "inference": ["predicted by MitoCurator intergenic-region scan"],
        },
    )
    record.features.append(feat)
    return {
        "seqid": record.id,
        "start": start_1b,
        "end": end_1b,
        "length": length,
        "at_content": round(at, 2),
        "feature_type": "misc_feature",
        "note": "putative AT-rich control region",
        "between": "intergenic_largest",
    }


def refine_annotation(config, input_gb, outdir):
    outdir = ensure_dir(outdir)
    record, _ = read_record(input_gb)

    summarize_expected_gene_set(record, outdir / "expected_gene_set.tsv")

    min_len = int(safe_get(config, ["refinement", "at_rich_min_len"], 500))
    min_at = float(safe_get(config, ["refinement", "at_rich_min_at"], 75.0))
    annotate_at = bool(safe_get(config, ["refinement", "annotate_at_rich"], True))

    added = add_at_rich_region(record, min_len=min_len, min_at=min_at) if annotate_at else None
    with open(outdir / "added_features.tsv", "w", encoding="utf-8") as out:
        out.write("seqid\tstart\tend\tlength\tat_content\tfeature_type\tnote\tbetween\n")
        if added:
            out.write("{seqid}\t{start}\t{end}\t{length}\t{at_content}\t{feature_type}\t{note}\t{between}\n".format(**added))

    out_gb = outdir / "refined.gb"
    write_record(record, out_gb, "genbank")
    return out_gb
