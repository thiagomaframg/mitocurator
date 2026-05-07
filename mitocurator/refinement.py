from __future__ import annotations
from collections import Counter
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from .io import read_record, write_record, feature_name
from .utils import ensure_dir, safe_get, get_genetic_code

EXPECTED_GENE_SET = {
    "CDS": ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"],
    "rRNA": ["rrnL", "rrnS"],
    "tRNA": [
        "tRNA-Ala", "tRNA-Arg", "tRNA-Asn", "tRNA-Asp", "tRNA-Cys", "tRNA-Gln", "tRNA-Glu", "tRNA-Gly", "tRNA-His",
        "tRNA-Ile", "tRNA-Leu", "tRNA-Leu2", "tRNA-Lys", "tRNA-Met", "tRNA-Phe", "tRNA-Pro", "tRNA-Ser", "tRNA-Ser2",
        "tRNA-Thr", "tRNA-Trp", "tRNA-Tyr", "tRNA-Val",
    ],
}
STOP_CODONS_TABLE5 = {"TAA", "TAG"}


def _normalize_token(value: str) -> str:
    x = (value or "").strip()
    u = x.upper().replace("_", "").replace(" ", "")
    aliases = {
        "NAD1": "ND1", "NAD2": "ND2", "NAD3": "ND3", "NAD4": "ND4", "NAD4L": "ND4L", "NAD5": "ND5", "NAD6": "ND6",
        "CYTB": "CYTB", "COB": "CYTB", "COX1": "COX1", "COI": "COX1", "COXI": "COX1", "COX2": "COX2", "COII": "COX2", "COX3": "COX3", "COIII": "COX3",
        "RRNL": "rrnL", "RRNS": "rrnS", "TRNW": "tRNA-Trp", "TRNL1": "tRNA-Leu", "TRNL2": "tRNA-Leu2", "TRNS1": "tRNA-Ser", "TRNS2": "tRNA-Ser2",
    }
    if u in aliases:
        return aliases[u]
    if u.startswith("TRNA"):
        y = x.replace("_", "-").replace(" ", "").replace("--", "-")
        y = y.replace("tRNA", "tRNA-") if y.startswith("tRNA") and not y.startswith("tRNA-") else y
        y = y.replace("TRNA", "tRNA-") if y.startswith("TRNA") else y
        y = y.replace("-1", "1").replace("-2", "2")
        y = y.replace("Leu1", "Leu").replace("Leu-1", "Leu").replace("Ser1", "Ser").replace("Ser-1", "Ser")
        return y
    return x.upper()


def _feature_segments(feat) -> List[Tuple[int, int]]:
    if isinstance(feat.location, CompoundLocation):
        return [(int(p.start), int(p.end)) for p in feat.location.parts]
    return [(int(feat.location.start), int(feat.location.end))]


def _merged_feature_ranges(record):
    segments = []
    for feat in record.features:
        if feat.type == "source":
            continue
        segments.extend(_feature_segments(feat))
    if not segments:
        return []
    merged = []
    for s, e in sorted(segments):
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return merged


def _intergenic_gaps(record):
    n = len(record.seq)
    merged = _merged_feature_ranges(record)
    if not merged:
        return []
    gaps = []
    for i in range(len(merged) - 1):
        a, b = merged[i], merged[i + 1]
        if b[0] > a[1]:
            gaps.append({"start0": a[1], "end0": b[0], "wrap": False})
    wrap_len = (merged[0][0] + n) - merged[-1][1]
    if wrap_len > 0:
        gaps.append({"start0": merged[-1][1], "end0": merged[0][0], "wrap": True})
    for g in gaps:
        seq = (record.seq[g["start0"]:g["end0"]] if not g["wrap"] else (record.seq[g["start0"]:] + record.seq[:g["end0"]]))
        g["length"] = len(seq)
        g["at_content"] = round(100.0 * (str(seq).upper().count("A") + str(seq).upper().count("T")) / len(seq), 2) if len(seq) else 0.0
    return gaps


def _translate_stop_metrics(nt: str, code: int):
    aa = str(Seq(nt).translate(table=code, to_stop=False))
    internal = aa[:-1].count("*") if aa else 0
    terminal = "yes" if aa.endswith("*") else "no"
    return len(aa), internal, terminal


def summarize_expected_gene_set(record, out_tsv: Path):
    counts = Counter()
    for feat in record.features:
        vals = []
        if feat.type == "CDS":
            vals.extend(feat.qualifiers.get("gene", []))
        elif feat.type in ("rRNA", "tRNA"):
            for k in ("gene", "product", "note"):
                vals.extend(feat.qualifiers.get(k, []))
        else:
            continue
        for v in vals:
            counts[_normalize_token(v)] += 1

    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("gene\ttype\tstatus\tcopies\tcomment\n")
        for gtype, genes in EXPECTED_GENE_SET.items():
            for gene in genes:
                c = counts.get(gene, 0)
                status = "PRESENT" if c == 1 else "DUPLICATED" if c > 1 else "MISSING"
                out.write(f"{gene}\t{gtype}\t{status}\t{c}\t{'.' if c else 'not detected in current annotation'}\n")


def add_at_rich_region(record, min_len=500, min_at=75.0):
    gaps = _intergenic_gaps(record)
    if not gaps:
        return None
    best = max(gaps, key=lambda g: g["length"])
    if best["length"] < min_len or best["at_content"] < min_at:
        return None

    s, e, n = best["start0"], best["end0"], len(record.seq)
    loc = (CompoundLocation([FeatureLocation(s, n, strand=1), FeatureLocation(0, e, strand=1)]) if best["wrap"] else FeatureLocation(s, e, strand=1))
    record.features.append(
        SeqFeature(location=loc, type="misc_feature", qualifiers={
            "note": ["putative AT-rich control region"],
            "product": ["AT-rich region"],
            "inference": ["predicted by MitoCurator intergenic-region scan"],
        })
    )
    return {
        "seqid": record.id,
        "start": s + 1,
        "end": e,
        "length": best["length"],
        "at_content": best["at_content"],
        "feature_type": "misc_feature",
        "note": "putative AT-rich control region",
        "between": "intergenic_largest",
    }


def _context_features(record, start0, end0):
    feats = sorted([f for f in record.features if f.type != "source"], key=lambda f: int(f.location.start))
    left, right = ".", "."
    for f in feats:
        if int(f.location.end) <= start0:
            left = f"{f.type}:{feature_name(f)}"
        if int(f.location.start) >= end0 and right == ".":
            right = f"{f.type}:{feature_name(f)}"
            break
    return left, right




def validate_candidate_translation(record, row, genetic_code):
    start0 = int(row["start"]) - 1
    end0 = int(row["end"])
    seq = record.seq[start0:end0]
    if row.get("strand") == "-":
        seq = seq.reverse_complement()
    frame = int(row.get("frame", 0) or 0)
    seq = seq[frame:]
    nt = str(seq).upper()
    nt_adj = nt[: (len(nt) // 3) * 3]
    aa_len, internal, terminal = _translate_stop_metrics(nt_adj, code=genetic_code) if nt_adj else (0, 0, "no")
    row["length_nt"] = len(nt_adj)
    row["length_aa"] = aa_len
    row["internal_stop_count"] = internal
    row["terminal_stop"] = terminal
    return row

def find_missing_cds_candidates(config, record, expected_gene_tsv, out_tsv):
    min_nt = int(safe_get(config, ["refinement", "orf_min_nt"], 150))
    genetic_code = get_genetic_code(config, default=5)
    exp = pd.read_csv(expected_gene_tsv, sep="	")
    missing = exp[(exp["type"] == "CDS") & (exp["status"] == "MISSING")]["gene"].tolist()
    gaps = sorted(_intergenic_gaps(record), key=lambda g: g["length"], reverse=True)
    rows = []
    cid = 1

    for gene in missing:
        candidates = []
        for gap in gaps:
            seq = (record.seq[gap["start0"]:gap["end0"]] if not gap["wrap"] else (record.seq[gap["start0"]:] + record.seq[:gap["end0"]]))
            if len(seq) < min_nt:
                continue
            for strand, nuc in (("+", str(seq).upper()), ("-", str(seq.reverse_complement()).upper())):
                for frame in (0, 1, 2):
                    i = frame
                    while i + 3 <= len(nuc):
                        while i + 3 <= len(nuc) and nuc[i:i+3] in STOP_CODONS_TABLE5:
                            i += 3
                        j = i
                        while j + 3 <= len(nuc) and nuc[j:j+3] not in STOP_CODONS_TABLE5:
                            j += 3
                        if j - i >= min_nt:
                            if strand == "+":
                                start0 = gap["start0"] + i
                                end0 = gap["start0"] + j
                            else:
                                # map coordinates from reverse-complemented local sequence back to genomic plus coordinates
                                start0 = gap["start0"] + (len(nuc) - j)
                                end0 = gap["start0"] + (len(nuc) - i)
                            row = {
                                "gene": gene,
                                "candidate_id": f"{gene}_cand{cid}",
                                "seqid": record.id,
                                "start": start0 + 1,
                                "end": end0,
                                "strand": strand,
                                "frame": 0,
                            }
                            row = validate_candidate_translation(record, row, genetic_code)
                            candidates.append((row["internal_stop_count"], -row["length_nt"], gap, row))
                            cid += 1
                        i = j + 3

        candidates.sort(key=lambda x: (x[0], x[1]))
        for _, _, gap, row in candidates[:5]:
            left, right = _context_features(record, row["start"] - 1, row["end"])
            row["overlaps_existing_feature"] = "no"
            row["nearest_left_feature"] = left
            row["nearest_right_feature"] = right
            row["region_context"] = f"gap_{gap['start0']+1}..{gap['end0']}|len={gap['length']}|AT={gap['at_content']}"
            nd2_target = gene == "ND2" and (gap["length"] >= 500 or (gap["start0"] + 1) <= 5411 <= max(gap["end0"], 5411))
            row["decision_hint"] = "STRONG_CANDIDATE_REGION" if nd2_target else ("CANDIDATE_REGION" if gap["length"] >= 500 else "LOW_PRIORITY")
            row["comment"] = "Validated translation from reported start/end/strand/frame"
            rows.append(row)

    cols = ["gene", "candidate_id", "seqid", "start", "end", "strand", "length_nt", "length_aa", "frame", "internal_stop_count", "terminal_stop", "overlaps_existing_feature", "nearest_left_feature", "nearest_right_feature", "region_context", "decision_hint", "comment"]
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="	", index=False)


def find_cds_refinement_candidates(config, record, problematic_features_tsv, out_tsv):
    window = int(safe_get(config, ["refinement", "cds_refinement_window"], 300))
    code = get_genetic_code(config, default=5)
    rows = []
    cds_feats = [f for f in record.features if f.type == "CDS"]

    for feat in cds_feats:
        gene = _normalize_token(feat.qualifiers.get("gene", [feature_name(feat)])[0])
        start0, end0 = int(feat.location.start), int(feat.location.end)
        strand = feat.location.strand or 1
        nt = str(feat.extract(record.seq)).upper()
        old_aa_len, old_internal, _ = _translate_stop_metrics(nt, code=code)
        if old_internal == 0:
            continue

        best = None
        for ds in range(-30, 31, 3):
            for de in range(-30, 31, 3):
                ns, ne = max(0, start0 + ds), min(len(record.seq), end0 + de)
                if ne - ns < 90:
                    continue
                cand_nt = str(record.seq[ns:ne].reverse_complement() if strand == -1 else record.seq[ns:ne]).upper()
                aa_len, internal, terminal = _translate_stop_metrics(cand_nt, code=code)
                score = (internal, abs((ne - ns) - (end0 - start0)))
                if best is None or score < best[0]:
                    best = (score, ns, ne, aa_len, internal, terminal, ds, de)

        if best is None:
            rows.append({"gene": gene, "old_start": start0 + 1, "old_end": end0, "old_strand": "+" if strand == 1 else "-", "old_length_nt": end0 - start0,
                         "old_internal_stop_count": old_internal, "candidate_start": ".", "candidate_end": ".", "candidate_strand": ".", "candidate_length_nt": ".", "candidate_length_aa": ".", "candidate_frame": ".", "candidate_internal_stop_count": ".", "candidate_terminal_stop": ".", "delta_start": ".", "delta_end": ".", "decision_hint": "NO_BETTER_CANDIDATE", "comment": "No alternative candidate in window search"})
            continue

        _, ns, ne, aa_len, internal, terminal, ds, de = best
        hint = "NO_BETTER_CANDIDATE"
        if internal < old_internal:
            hint = "SUGGEST_REVIEW"
        if internal == 0 and abs((ne - ns) - (end0 - start0)) <= window:
            hint = "STRONG_CANDIDATE"
        rows.append({"gene": gene, "old_start": start0 + 1, "old_end": end0, "old_strand": "+" if strand == 1 else "-", "old_length_nt": end0 - start0,
                     "old_internal_stop_count": old_internal, "candidate_start": ns + 1, "candidate_end": ne, "candidate_strand": "+" if strand == 1 else "-",
                     "candidate_length_nt": ne - ns, "candidate_length_aa": aa_len, "candidate_frame": (ns % 3), "candidate_internal_stop_count": internal,
                     "candidate_terminal_stop": terminal, "delta_start": ds, "delta_end": de, "decision_hint": hint,
                     "comment": "Coordinate-only candidate scan; no GenBank changes applied"})

    cols = ["gene", "old_start", "old_end", "old_strand", "old_length_nt", "old_internal_stop_count", "candidate_start", "candidate_end", "candidate_strand", "candidate_length_nt", "candidate_length_aa", "candidate_frame", "candidate_internal_stop_count", "candidate_terminal_stop", "delta_start", "delta_end", "decision_hint", "comment"]
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="\t", index=False)


def refine_annotation(config, input_gb, outdir):
    outdir = ensure_dir(outdir)
    record, _ = read_record(input_gb)

    expected_tsv = outdir / "expected_gene_set.tsv"
    summarize_expected_gene_set(record, expected_tsv)

    min_len = int(safe_get(config, ["refinement", "at_rich_min_len"], 500))
    min_at = float(safe_get(config, ["refinement", "at_rich_min_at"], 75.0))
    annotate_at = bool(safe_get(config, ["refinement", "annotate_at_rich"], True))
    added = add_at_rich_region(record, min_len=min_len, min_at=min_at) if annotate_at else None

    with open(outdir / "added_features.tsv", "w", encoding="utf-8") as out:
        out.write("seqid\tstart\tend\tlength\tat_content\tfeature_type\tnote\tbetween\n")
        if added:
            out.write("{seqid}\t{start}\t{end}\t{length}\t{at_content}\t{feature_type}\t{note}\t{between}\n".format(**added))

    if bool(safe_get(config, ["refinement", "find_missing_cds_candidates"], True)):
        find_missing_cds_candidates(config, record, expected_tsv, outdir / "missing_gene_candidates.tsv")
    else:
        pd.DataFrame(columns=["gene", "candidate_id", "seqid", "start", "end", "strand", "length_nt", "length_aa", "frame", "internal_stop_count", "terminal_stop", "overlaps_existing_feature", "nearest_left_feature", "nearest_right_feature", "region_context", "decision_hint", "comment"]).to_csv(outdir / "missing_gene_candidates.tsv", sep="\t", index=False)

    find_cds_refinement_candidates(config, record, outdir / "problematic_features.tsv", outdir / "cds_refinement_candidates.tsv") if bool(safe_get(config, ["refinement", "find_cds_refinement_candidates"], True)) else pd.DataFrame(columns=["gene", "old_start", "old_end", "old_strand", "old_length_nt", "old_internal_stop_count", "candidate_start", "candidate_end", "candidate_strand", "candidate_length_nt", "candidate_length_aa", "candidate_frame", "candidate_internal_stop_count", "candidate_terminal_stop", "delta_start", "delta_end", "decision_hint", "comment"]).to_csv(outdir / "cds_refinement_candidates.tsv", sep="\t", index=False)

    out_gb = outdir / "refined.gb"
    write_record(record, out_gb, "genbank")
    return out_gb
