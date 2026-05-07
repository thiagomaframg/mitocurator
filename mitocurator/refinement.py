from __future__ import annotations
from collections import Counter
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
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

DEFAULT_METAZOA = {
    "CDS": ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"],
    "rRNA": ["rrnL", "rrnS"],
    "tRNA": [
        "tRNA-Ala", "tRNA-Arg", "tRNA-Asn", "tRNA-Asp", "tRNA-Cys", "tRNA-Gln", "tRNA-Glu", "tRNA-Gly", "tRNA-His", "tRNA-Ile",
        "tRNA-Leu", "tRNA-Leu2", "tRNA-Lys", "tRNA-Met", "tRNA-Phe", "tRNA-Pro", "tRNA-Ser", "tRNA-Ser2", "tRNA-Thr", "tRNA-Trp", "tRNA-Tyr", "tRNA-Val",
    ],
}


def _get_expected_gene_cfg(config):
    raw = safe_get(config, ["refinement", "expected_gene_set"], "insect_mito")
    if isinstance(raw, str):
        return {"profile": raw, "custom_file": None}
    if isinstance(raw, dict):
        return {"profile": raw.get("profile", "insect_mito"), "custom_file": raw.get("custom_file")}
    return {"profile": "insect_mito", "custom_file": None}


def _load_expected_profile(config):
    cfg = _get_expected_gene_cfg(config)
    profile = str(cfg["profile"] or "insect_mito")
    if profile in ("metazoa_mito", "insect_mito", "vertebrate_mito"):
        rows = []
        for t, genes in DEFAULT_METAZOA.items():
            for g in genes:
                rows.append({"gene": g, "type": t, "required": True, "aliases": ""})
        return profile, rows
    if profile == "minimal_mito":
        return profile, []
    if profile == "custom":
        cfile = cfg.get("custom_file")
        if not cfile:
            raise RuntimeError("expected_gene_set.profile=custom requires expected_gene_set.custom_file")
        df = pd.read_csv(cfile, sep="	")
        req_cols = {"gene", "type", "required", "aliases"}
        if not req_cols.issubset(df.columns):
            raise RuntimeError("custom expected gene file must include: gene, type, required, aliases")
        rows = []
        for _, r in df.iterrows():
            rows.append({"gene": str(r["gene"]), "type": str(r["type"]), "required": str(r["required"]).lower() in ("1", "true", "yes"), "aliases": str(r.get("aliases", "") or "")})
        return profile, rows
    return "insect_mito", _load_expected_profile({"refinement": {"expected_gene_set": "insect_mito"}})[1]



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


def summarize_expected_gene_set(config, record, out_tsv: Path):
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

    profile, expected_rows = _load_expected_profile(config)
    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("gene\ttype\tstatus\tcopies\tprofile\trequired\tcomment\n")
        if profile == "minimal_mito":
            for g, c in sorted(counts.items()):
                out.write(f"{g}\tobserved\tPRESENT\t{c}\t{profile}\tfalse\tminimal profile: observed-only summary\n")
            return
        for row in expected_rows:
            gene = _normalize_token(row["gene"])
            c = counts.get(gene, 0)
            required = bool(row.get("required", True))
            if not required:
                status = "OPTIONAL"
            else:
                status = "PRESENT" if c == 1 else "DUPLICATED" if c > 1 else "MISSING"
            comment = "." if c else ("optional gene not detected" if not required else "not detected in current annotation")
            out.write(f"{gene}\t{row['type']}\t{status}\t{c}\t{profile}\t{str(required).lower()}\t{comment}\n")


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




def _reference_path_from_config(config):
    return (
        safe_get(config, ["reference", "genbank"], None)
        or safe_get(config, ["mitofinder", "reference_gb"], None)
        or safe_get(config, ["mitofinder", "reference"], None)
    )


def _extract_reference_cds_proteins(config):
    ref_path = _reference_path_from_config(config)
    if not ref_path or not Path(ref_path).exists():
        return None, {}
    ref_record, _ = read_record(ref_path)
    code = get_genetic_code(config, default=5)
    prots = {}
    for feat in ref_record.features:
        if feat.type != "CDS":
            continue
        gene = _normalize_token(feat.qualifiers.get("gene", [feature_name(feat)])[0])
        nt = str(feat.extract(ref_record.seq)).upper()
        nt = nt[: (len(nt) // 3) * 3]
        if not nt:
            continue
        aa = str(Seq(nt).translate(table=code, to_stop=False)).replace("*", "")
        if aa and (gene not in prots or len(aa) > len(prots[gene])):
            prots[gene] = aa
    return ref_path, prots


def _pairwise_metrics(query_aa, ref_aa):
    if not query_aa or not ref_aa:
        return 0.0, 0.0, 0.0, 0.0
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    al = aligner.align(query_aa, ref_aa)[0]
    score = float(al.score)
    q_aln, r_aln = al[0], al[1]

    overlap = sum(1 for a, b in zip(q_aln, r_aln) if a != "-" and b != "-")
    matches = sum(1 for a, b in zip(q_aln, r_aln) if a == b and a != "-" and b != "-")
    pid = (100.0 * matches / overlap) if overlap else 0.0

    # Coverage based on truly aligned non-gap overlap only.
    cov_q = (100.0 * overlap / len(query_aa)) if query_aa else 0.0
    cov_r = (100.0 * overlap / len(ref_aa)) if ref_aa else 0.0

    # Approximate robust cap to avoid artificial full coverage from global alignment padding.
    approx_ref = (100.0 * min(len(query_aa), len(ref_aa)) / len(ref_aa)) if ref_aa else 0.0
    approx_q = (100.0 * min(len(query_aa), len(ref_aa)) / len(query_aa)) if query_aa else 0.0
    cov_q = min(cov_q, approx_q)
    cov_r = min(cov_r, approx_ref)

    return round(score, 2), round(pid, 2), round(cov_q, 2), round(cov_r, 2)


def generate_reference_similarity_candidates(config, record, missing_candidates_tsv, out_tsv):
    cols = ["gene","candidate_id","seqid","start","end","strand","length_nt","candidate_aa_length","reference_gene","reference_aa_length","length_ratio","percent_identity","aligned_coverage_candidate","aligned_coverage_reference","alignment_score","internal_stop_count","terminal_stop","region_context","decision_hint","comment"]
    ref_path, ref_prot = _extract_reference_cds_proteins(config)
    if not Path(missing_candidates_tsv).exists():
        pd.DataFrame(columns=cols).to_csv(out_tsv, sep="	", index=False)
        return
    cand_df = pd.read_csv(missing_candidates_tsv, sep="	")
    if cand_df.empty:
        pd.DataFrame(columns=cols).to_csv(out_tsv, sep="	", index=False)
        return
    code = get_genetic_code(config, default=5)
    min_id = float(safe_get(config, ["refinement", "reference_similarity_min_identity"], 25.0))
    rows = []
    for _, r in cand_df.iterrows():
        row = {k: r.get(k, ".") for k in ["gene","candidate_id","seqid","start","end","strand","length_nt","internal_stop_count","terminal_stop","region_context"]}
        gene = _normalize_token(str(r["gene"]))
        start0, end0 = int(r["start"])-1, int(r["end"])
        seq = record.seq[start0:end0]
        if str(r["strand"]) == "-":
            seq = seq.reverse_complement()
        frame = int(r.get("frame", 0) or 0)
        seq = seq[frame:]
        nt = str(seq).upper()
        nt = nt[: (len(nt)//3)*3]
        aa = str(Seq(nt).translate(table=code, to_stop=False)) if nt else ""
        aa_nostop = aa.replace("*", "")
        row["candidate_aa_length"] = len(aa_nostop)
        ref_aa = ref_prot.get(gene)
        if not ref_path:
            row.update({"reference_gene": ".", "reference_aa_length": ".", "length_ratio": ".", "percent_identity": ".", "aligned_coverage_candidate": ".", "aligned_coverage_reference": ".", "alignment_score": ".", "decision_hint": "NO_REFERENCE_GENE", "comment": "reference_gb not provided"})
        elif ref_aa is None:
            row.update({"reference_gene": gene, "reference_aa_length": ".", "length_ratio": ".", "percent_identity": ".", "aligned_coverage_candidate": ".", "aligned_coverage_reference": ".", "alignment_score": ".", "decision_hint": "NO_REFERENCE_GENE", "comment": "gene not present in reference CDS"})
        else:
            score, pid, covq, covr = _pairwise_metrics(aa_nostop, ref_aa)
            ratio = round((len(aa_nostop)/len(ref_aa)), 3) if len(ref_aa) else 0
            hint = "WEAK_REF_MATCH"
            if int(r.get("internal_stop_count", 0)) > 0:
                hint = "HAS_INTERNAL_STOP"
            elif pid >= 40 and covr >= 70 and ratio >= 0.70:
                hint = "STRONG_REF_MATCH"
            elif pid >= min_id and covr >= 30:
                hint = "PARTIAL_REF_MATCH"
            row.update({"reference_gene": gene, "reference_aa_length": len(ref_aa), "length_ratio": ratio, "percent_identity": pid, "aligned_coverage_candidate": covq, "aligned_coverage_reference": covr, "alignment_score": score, "decision_hint": hint, "comment": "protein alignment to reference CDS; overlap-based coverage (non-gap only)"})
        rows.append(row)
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="	", index=False)


def generate_problematic_cds_reference_check(config, record, out_tsv):
    cols = ["gene","seqid","start","end","strand","length_nt","candidate_aa_length","reference_gene","reference_aa_length","length_ratio","percent_identity","aligned_coverage_candidate","aligned_coverage_reference","alignment_score","internal_stop_count","terminal_stop","decision_hint","comment"]
    ref_path, ref_prot = _extract_reference_cds_proteins(config)
    code = get_genetic_code(config, default=5)
    rows = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        gene = _normalize_token(feat.qualifiers.get("gene", [feature_name(feat)])[0])
        start0, end0 = int(feat.location.start), int(feat.location.end)
        nt = str(feat.extract(record.seq)).upper()
        nt3 = nt[: (len(nt)//3)*3]
        aa_len, internal, terminal = _translate_stop_metrics(nt3, code=code) if nt3 else (0,0,"no")
        if internal <= 0:
            continue
        aa = str(Seq(nt3).translate(table=code, to_stop=False)).replace("*", "") if nt3 else ""
        base = {"gene": gene, "seqid": record.id, "start": start0+1, "end": end0, "strand": "+" if feat.location.strand != -1 else "-", "length_nt": len(nt3), "candidate_aa_length": len(aa), "internal_stop_count": internal, "terminal_stop": terminal}
        if not ref_path:
            base.update({"reference_gene": ".", "reference_aa_length": ".", "length_ratio": ".", "percent_identity": ".", "aligned_coverage_candidate": ".", "aligned_coverage_reference": ".", "alignment_score": ".", "decision_hint": "NO_REFERENCE_GENE", "comment": "reference_gb not provided"})
        elif gene not in ref_prot:
            base.update({"reference_gene": gene, "reference_aa_length": ".", "length_ratio": ".", "percent_identity": ".", "aligned_coverage_candidate": ".", "aligned_coverage_reference": ".", "alignment_score": ".", "decision_hint": "NO_REFERENCE_GENE", "comment": "gene not present in reference CDS"})
        else:
            refaa = ref_prot[gene]
            score, pid, covq, covr = _pairwise_metrics(aa, refaa)
            ratio = round((len(aa)/len(refaa)), 3) if len(refaa) else 0
            if pid >= 40 and covr >= 50:
                hint = "POSSIBLE_SEQUENCE_ERROR"
            elif covr >= 30:
                hint = "PARTIAL_MATCH_WITH_STOPS"
            else:
                hint = "LIKELY_COORDINATE_OR_FRAME_PROBLEM"
            base.update({"reference_gene": gene, "reference_aa_length": len(refaa), "length_ratio": ratio, "percent_identity": pid, "aligned_coverage_candidate": covq, "aligned_coverage_reference": covr, "alignment_score": score, "decision_hint": hint, "comment": "problematic CDS translated and compared to reference"})
        rows.append(base)
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="	", index=False)



def _problematic_cds_features(record, genetic_code):
    rows = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        gene = _normalize_token(feat.qualifiers.get("gene", [feature_name(feat)])[0])
        nt = str(feat.extract(record.seq)).upper()
        nt3 = nt[: (len(nt)//3)*3]
        aa = str(Seq(nt3).translate(table=genetic_code, to_stop=False)) if nt3 else ""
        internal_positions = [i+1 for i, r in enumerate(aa[:-1]) if r == "*"]
        if not internal_positions:
            continue
        rows.append((feat, gene, nt3, aa, internal_positions))
    return rows


def generate_problematic_cds_stop_context(config, record, out_tsv):
    cols = ["gene","seqid","cds_start","cds_end","strand","length_nt","length_aa","stop_aa_position","genomic_codon_start","genomic_codon_end","codon","translation_table","previous_aa","next_aa","comment"]
    code = get_genetic_code(config, default=5)
    rows = []
    for feat, gene, nt3, aa, positions in _problematic_cds_features(record, code):
        start0, end0 = int(feat.location.start), int(feat.location.end)
        strand = feat.location.strand or 1
        for pos in positions:
            codon_idx = (pos - 1) * 3
            codon = nt3[codon_idx:codon_idx+3]
            if strand == 1:
                gstart = start0 + codon_idx + 1
                gend = gstart + 2
            else:
                gend = end0 - codon_idx
                gstart = gend - 2
            rows.append({
                "gene": gene, "seqid": record.id, "cds_start": start0 + 1, "cds_end": end0,
                "strand": "+" if strand == 1 else "-", "length_nt": len(nt3), "length_aa": len(aa),
                "stop_aa_position": pos, "genomic_codon_start": gstart, "genomic_codon_end": gend,
                "codon": codon, "translation_table": code,
                "previous_aa": aa[pos-2] if pos > 1 else ".", "next_aa": aa[pos] if pos < len(aa) else ".",
                "comment": "internal stop codon context",
            })
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="	", index=False)


def _map_candidate_to_reference_positions(candidate_aa, reference_aa):
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    aln = aligner.align(candidate_aa, reference_aa)[0]
    ca, ra = aln[0], aln[1]
    cmap = {}
    cpos = rpos = 0
    for c, r in zip(ca, ra):
        if c != "-":
            cpos += 1
        if r != "-":
            rpos += 1
        if c != "-":
            cmap[cpos] = rpos if r != "-" else None
    return cmap


def generate_problematic_cds_reference_alignment(config, record, out_tsv):
    cols = ["gene","candidate_aa_position","candidate_residue","reference_gene","reference_aa_position","reference_residue","candidate_context","reference_context","comment"]
    code = get_genetic_code(config, default=5)
    _, ref_prot = _extract_reference_cds_proteins(config)
    rows = []
    for _, gene, _, aa, positions in _problematic_cds_features(record, code):
        refaa = ref_prot.get(gene)
        if not refaa:
            for pos in positions:
                rows.append({"gene": gene, "candidate_aa_position": pos, "candidate_residue": "*", "reference_gene": gene, "reference_aa_position": ".", "reference_residue": ".", "candidate_context": aa[max(0,pos-4):pos+3], "reference_context": ".", "comment": "reference gene not found"})
            continue
        mp = _map_candidate_to_reference_positions(aa, refaa)
        for pos in positions:
            rpos = mp.get(pos)
            rres = refaa[rpos-1] if rpos and rpos <= len(refaa) else "."
            rctx = refaa[max(0,(rpos or 1)-4):(rpos or 1)+3] if rpos else "."
            rows.append({"gene": gene, "candidate_aa_position": pos, "candidate_residue": "*", "reference_gene": gene, "reference_aa_position": rpos if rpos else ".", "reference_residue": rres, "candidate_context": aa[max(0,pos-4):pos+3], "reference_context": rctx, "comment": "mapped by global protein alignment"})
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="	", index=False)


def write_missing_gene_candidate_proteins(config, record, candidates_tsv, out_faa):
    code = get_genetic_code(config, default=5)
    lines = []
    if Path(candidates_tsv).exists():
        df = pd.read_csv(candidates_tsv, sep="	")
        for _, r in df.iterrows():
            s0, e0 = int(r["start"])-1, int(r["end"])
            seq = record.seq[s0:e0]
            if str(r["strand"]) == "-":
                seq = seq.reverse_complement()
            frame = int(r.get("frame",0) or 0)
            seq = seq[frame:]
            nt = str(seq).upper()
            nt3 = nt[:(len(nt)//3)*3]
            aa = str(Seq(nt3).translate(table=code, to_stop=False)).replace("*", "") if nt3 else ""
            header = f">{r['gene']}|{r['candidate_id']}|{int(r['start'])}..{int(r['end'])}|{r['strand']}|length_aa={len(aa)}|{r.get('decision_hint','.') }"
            lines.extend([header, aa])

    Path(out_faa).write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")

def write_problematic_cds_proteins(config, record, out_faa):
    code = get_genetic_code(config, default=5)
    lines = []
    for feat, gene, nt3, aa, positions in _problematic_cds_features(record, code):
        s, e = int(feat.location.start)+1, int(feat.location.end)
        strand = "+" if (feat.location.strand or 1) == 1 else "-"
        header = f">{gene}|{s}..{e}|{strand}|length_aa={len(aa)}|internal_stops={len(positions)}"
        lines.extend([header, aa])

    Path(out_faa).write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def _read_tsv_if_exists(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(p, sep="	")


def _hint_rank(h):
    order = {"STRONG_REF_MATCH": 3, "PARTIAL_REF_MATCH": 2, "WEAK_REF_MATCH": 1}
    return order.get(str(h), 0)


def generate_curation_recommendations(config, outdir):
    outdir = Path(outdir)
    expected = _read_tsv_if_exists(outdir / "expected_gene_set.tsv")
    refcand = _read_tsv_if_exists(outdir / "reference_similarity_candidates.tsv")
    prob = _read_tsv_if_exists(outdir / "problematic_cds_reference_check.tsv")
    stops = _read_tsv_if_exists(outdir / "problematic_cds_stop_context.tsv")

    cols = ["target","target_type","issue_type","status","best_candidate","best_candidate_coordinates","best_candidate_strand","best_candidate_length_aa","reference_gene","reference_length_aa","length_ratio","percent_identity","aligned_coverage_reference","internal_stop_count","stop_count","stop_positions","evidence_summary","recommendation","priority","comment"]
    rows = []

    # Missing genes
    if not expected.empty:
        miss = expected[(expected.get("type") == "CDS") & (expected.get("status") == "MISSING")]
        for _, m in miss.iterrows():
            gene = str(m["gene"])
            cands = refcand[refcand["gene"] == gene] if not refcand.empty and "gene" in refcand.columns else pd.DataFrame()
            if cands.empty:
                rows.append({"target": gene, "target_type": "CDS", "issue_type": "MISSING_GENE", "status": "MISSING", "best_candidate": ".", "best_candidate_coordinates": ".", "best_candidate_strand": ".", "best_candidate_length_aa": ".", "reference_gene": gene, "reference_length_aa": ".", "length_ratio": ".", "percent_identity": ".", "aligned_coverage_reference": ".", "internal_stop_count": ".", "stop_count": 0, "stop_positions": ".", "evidence_summary": "No candidate found in reference similarity table", "recommendation": "NO_CANDIDATE_FOUND", "priority": "LOW", "comment": "diagnostic-only"})
                continue
            cands = cands.copy()
            for c in ["aligned_coverage_reference","length_ratio","percent_identity"]:
                cands[c] = pd.to_numeric(cands[c], errors="coerce").fillna(0)
            cands["_rank"] = cands["decision_hint"].map(_hint_rank).fillna(0)
            best = cands.sort_values(["_rank","aligned_coverage_reference","length_ratio","percent_identity"], ascending=False).iloc[0]
            rec = "NO_CANDIDATE_FOUND"
            pri = "LOW"
            if best["decision_hint"] == "PARTIAL_REF_MATCH":
                rec, pri = "REVIEW_PARTIAL_CANDIDATE", "MEDIUM"
            elif best["decision_hint"] == "STRONG_REF_MATCH":
                rec, pri = "REVIEW_FOR_MANUAL_ANNOTATION", "HIGH"
            rows.append({"target": gene, "target_type": "CDS", "issue_type": "MISSING_GENE", "status": "MISSING", "best_candidate": best.get("candidate_id","."), "best_candidate_coordinates": f"{best.get('start','.') }..{best.get('end','.')}", "best_candidate_strand": best.get("strand","."), "best_candidate_length_aa": best.get("candidate_aa_length","."), "reference_gene": best.get("reference_gene",gene), "reference_length_aa": best.get("reference_aa_length","."), "length_ratio": best.get("length_ratio","."), "percent_identity": best.get("percent_identity","."), "aligned_coverage_reference": best.get("aligned_coverage_reference","."), "internal_stop_count": best.get("internal_stop_count","."), "stop_count": int(best.get("internal_stop_count",0) or 0), "stop_positions": ".", "evidence_summary": f"best_hint={best.get('decision_hint','.')}; cov_ref={best.get('aligned_coverage_reference','.')}; ratio={best.get('length_ratio','.')}", "recommendation": rec, "priority": pri, "comment": "manual curation only; no automatic annotation"})

    # Problematic CDS with internal stops
    if not prob.empty:
        for _, r in prob.iterrows():
            gene = str(r["gene"])
            rratio = float(r.get("length_ratio", 0) or 0)
            pid = float(r.get("percent_identity", 0) or 0)
            covr = float(r.get("aligned_coverage_reference", 0) or 0)
            istops = int(r.get("internal_stop_count", 0) or 0)
            srows = stops[stops["gene"] == gene] if not stops.empty and "gene" in stops.columns else pd.DataFrame()
            stop_pos = sorted(srows["stop_aa_position"].tolist()) if not srows.empty else []
            clustered = (len(stop_pos) >= 3 and max(stop_pos) - min(stop_pos) <= 30) if stop_pos else False
            if rratio >= 0.90 and covr >= 80 and istops > 0:
                rec, pri = "CHECK_READ_SUPPORT_FOR_LOCAL_ERROR", "HIGH"
            elif istops >= 3 and clustered:
                rec, pri = "CHECK_FRAMESHIFT_OR_LOCAL_INDEL", "HIGH"
            elif rratio < 0.70 and pid >= 50:
                rec, pri = "SEARCH_EXTENDED_REGION_OR_REANNOTATE_BOUNDARIES", "MEDIUM"
            elif pid < 40 or covr < 30:
                rec, pri = "POSSIBLE_MISANNOTATION", "LOW"
            else:
                rec, pri = "MANUAL_REVIEW", "MEDIUM"
            rows.append({"target": gene, "target_type": "CDS", "issue_type": "INTERNAL_STOP", "status": "PROBLEMATIC", "best_candidate": ".", "best_candidate_coordinates": f"{r.get('start','.') }..{r.get('end','.')}", "best_candidate_strand": r.get("strand","."), "best_candidate_length_aa": r.get("candidate_aa_length","."), "reference_gene": r.get("reference_gene",gene), "reference_length_aa": r.get("reference_aa_length","."), "length_ratio": r.get("length_ratio","."), "percent_identity": r.get("percent_identity","."), "aligned_coverage_reference": r.get("aligned_coverage_reference","."), "internal_stop_count": istops, "stop_count": len(stop_pos) if stop_pos else istops, "stop_positions": ",".join(map(str, stop_pos)) if stop_pos else ".", "evidence_summary": f"pid={pid}; cov_ref={covr}; ratio={rratio}; internal_stops={istops}", "recommendation": rec, "priority": pri, "comment": "manual curation only; no automatic CDS fix"})

    rec_df = pd.DataFrame(rows, columns=cols)
    rec_df.to_csv(outdir / "curation_recommendations.tsv", sep="	", index=False)

    high = rec_df[rec_df["priority"] == "HIGH"] if not rec_df.empty else pd.DataFrame()
    missing_n = int((rec_df["issue_type"] == "MISSING_GENE").sum()) if not rec_df.empty else 0
    stops_n = int((rec_df["issue_type"] == "INTERNAL_STOP").sum()) if not rec_df.empty else 0
    with open(outdir / "curation_recommendations.md", "w", encoding="utf-8") as md:
        md.write("# MitoCurator curation recommendations\n\n")
        md.write(f"- Total recommendations: {len(rec_df)}\n")
        md.write(f"- Missing-gene targets: {missing_n}\n")
        md.write(f"- Problematic CDS targets: {stops_n}\n\n")
        md.write("## Missing genes\n\n")
        if rec_df.empty or missing_n == 0:
            md.write("No missing-gene recommendations.\n")
        else:
            for _, r in rec_df[rec_df["issue_type"] == "MISSING_GENE"].iterrows():
                md.write(f"- **{r['target']}**: {r['recommendation']} (priority: {r['priority']}) — {r['evidence_summary']}\n")
        md.write("\n## CDS with internal stops\n\n")
        if rec_df.empty or stops_n == 0:
            md.write("No internal-stop CDS recommendations.\n")
        else:
            for _, r in rec_df[rec_df["issue_type"] == "INTERNAL_STOP"].iterrows():
                md.write(f"- **{r['target']}**: {r['recommendation']} (priority: {r['priority']}) — {r['evidence_summary']}\n")
        md.write("\n## Prioritized recommendations (HIGH)\n\n")
        if high.empty:
            md.write("No HIGH-priority recommendations.\n")
        else:
            for _, r in high.iterrows():
                md.write(f"- {r['target']} ({r['issue_type']}): {r['recommendation']}\n")
        md.write("\n> No automatic GenBank correction was applied; this report is diagnostic-only.\n")


def refine_annotation(config, input_gb, outdir):
    outdir = ensure_dir(outdir)
    record, _ = read_record(input_gb)

    expected_tsv = outdir / "expected_gene_set.tsv"
    summarize_expected_gene_set(config, record, expected_tsv)

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

    generate_problematic_cds_stop_context(config, record, outdir / "problematic_cds_stop_context.tsv")
    generate_problematic_cds_reference_alignment(config, record, outdir / "problematic_cds_reference_alignment.tsv")

    if bool(safe_get(config, ["refinement", "compare_candidates_to_reference"], False)):
        generate_reference_similarity_candidates(
            config,
            record,
            outdir / "missing_gene_candidates.tsv",
            outdir / "reference_similarity_candidates.tsv",
        )

    if bool(safe_get(config, ["refinement", "compare_problematic_cds_to_reference"], False)):
        generate_problematic_cds_reference_check(
            config,
            record,
            outdir / "problematic_cds_reference_check.tsv",
        )


    if bool(safe_get(config, ["refinement", "compare_candidates_to_reference"], False)):
        generate_reference_similarity_candidates(
            config,
            record,
            outdir / "missing_gene_candidates.tsv",
            outdir / "reference_similarity_candidates.tsv",
        )

    if bool(safe_get(config, ["refinement", "compare_problematic_cds_to_reference"], False)):
        generate_problematic_cds_reference_check(
            config,
            record,
            outdir / "problematic_cds_reference_check.tsv",
        )

    write_missing_gene_candidate_proteins(config, record, outdir / "missing_gene_candidates.tsv", outdir / "missing_gene_candidate_proteins.faa")
    write_problematic_cds_proteins(config, record, outdir / "problematic_cds_proteins.faa")
    generate_curation_recommendations(config, outdir)

    out_gb = outdir / "refined.gb"
    write_record(record, out_gb, "genbank")
    return out_gb
