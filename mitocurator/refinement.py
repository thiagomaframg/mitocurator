from __future__ import annotations
from collections import Counter
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from .io import read_record, write_record, feature_name
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


def _translate_stop_metrics(nt: str, code: int = 5):
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


def _normalize_trna_name(name: str) -> str:
    return name.lower().replace("trna", "").replace("-", "").replace("_", "").strip()


def _find_trna(record, name: str):
    """Find a tRNA feature matching name; case-insensitive, ignores 'tRNA' prefix and separators."""
    norm = _normalize_trna_name(name)
    for f in record.features:
        if f.type != "tRNA":
            continue
        for qual in ("gene", "product", "note"):
            for v in f.qualifiers.get(qual, []):
                if _normalize_trna_name(v) == norm:
                    return f
    return None


def _gaps_adjacent_to(gaps, pos, tolerance=10):
    """Return all gaps whose start0 or end0 is within tolerance bp of pos."""
    return [g for g in gaps
            if abs(g["start0"] - pos) <= tolerance or abs(g["end0"] - pos) <= tolerance]


def _find_control_region_by_flanks(record, flanks, gaps):
    """Locate the control region gap using tRNA flanking markers.

    flanks is a list of [five_prime_name, three_prime_name] pairs, tried in order.
    For Hymenoptera the canonical pair is ["tRNA-Ile", "tRNA-Gln"].
    For vertebrates the typical pair is ["tRNA-Pro", "tRNA-Phe"] (D-loop context).
    Other metazoan groups may use different pairs — check the relevant literature.

    Strategy per pair:
      - Both found: control region is the gap cleanly between them.
      - Only one found: largest gap adjacent (≤10 bp) to that tRNA.
      - Neither found: try next pair; return (None, None) if all pairs exhausted.

    Returns (gap_dict, method_note) or (None, None).
    """
    for five_prime, three_prime in flanks:
        t5 = _find_trna(record, five_prime)
        t3 = _find_trna(record, three_prime)

        if t5 is not None and t3 is not None:
            e5 = int(t5.location.end)
            s3 = int(t3.location.start)
            e3 = int(t3.location.end)
            s5 = int(t5.location.start)
            for g in gaps:
                if ((abs(g["start0"] - e5) <= 10 and abs(g["end0"] - s3) <= 10) or
                        (abs(g["start0"] - e3) <= 10 and abs(g["end0"] - s5) <= 10)):
                    return g, f"tRNA flanks: {five_prime} + {three_prime}"
            # Flanks found but no single gap cleanly between them; fall to single-anchor.

        if t3 is not None:
            adj = _gaps_adjacent_to(gaps, int(t3.location.end)) + \
                  _gaps_adjacent_to(gaps, int(t3.location.start))
            adj = list({id(g): g for g in adj}.values())  # deduplicate
            if adj:
                # Selects the largest adjacent gap. Validated for M. capixaba (HiFi, Hymenoptera):
                # control region (2136 bp) was unambiguously larger than the only other adjacent
                # gap (29 bp, between rrnS and tRNA-Gln). If a future case presents multiple gaps
                # of similar size near a single anchor, this choice becomes arbitrary — should be
                # treated as ambiguous (confidence=low or rejection). Known limitation; not
                # implemented, flagged for revision. See docs/mitocurator_dev_brief.md §Limitações.
                return max(adj, key=lambda g: g["length"]), \
                       f"tRNA flank: {three_prime} (single anchor; {five_prime} not found)"

        if t5 is not None:
            adj = _gaps_adjacent_to(gaps, int(t5.location.end)) + \
                  _gaps_adjacent_to(gaps, int(t5.location.start))
            adj = list({id(g): g for g in adj}.values())
            if adj:
                # Same "largest adjacent gap" heuristic — see comment above.
                return max(adj, key=lambda g: g["length"]), \
                       f"tRNA flank: {five_prime} (single anchor; {three_prime} not found)"

    return None, None


def _tandem_repeat_score(seq: str, block_size: int = 55, min_identity: float = 0.75,
                         min_copies: int = 3) -> tuple[bool, dict]:
    """Detect tandem repeats by pairwise comparison of non-overlapping blocks.

    Splits seq into non-overlapping blocks of block_size bp, then for each block
    counts how many other blocks have pairwise identity >= min_identity (including
    itself, so the minimum count for a match is 1). Returns True if any block
    finds >= min_copies similar blocks, i.e. the repeat unit appears >= min_copies
    times in the candidate window.
    """
    seq = seq.upper()
    n = len(seq)
    if n < block_size * min_copies:
        return False, {"max_copies_found": 0}
    blocks = [seq[i:i + block_size] for i in range(0, n, block_size)
              if i + block_size <= n]
    if len(blocks) < min_copies:
        return False, {"max_copies_found": 0}
    max_copies = 0
    for query in blocks:
        copies = sum(
            sum(a == b for a, b in zip(query, blk)) / block_size >= min_identity
            for blk in blocks
        )
        if copies > max_copies:
            max_copies = copies
    return max_copies >= min_copies, {"max_copies_found": max_copies}


def add_at_rich_region(record, min_len=500, min_at=75.0,
                       control_region_flanks=None,
                       repeat_block_size=55, repeat_min_identity=0.75, repeat_min_copies=3):
    """Annotate the most likely mitochondrial control region as a misc_feature.

    Detection uses three signals in priority order:

    1. tRNA-flank detection (primary, biological prior):
       control_region_flanks is a list of [five_prime_tRNA, three_prime_tRNA] name pairs
       tried in order. The default ["tRNA-Ile", "tRNA-Gln"] is appropriate for Hymenoptera.
       Vertebrates typically use ["tRNA-Pro", "tRNA-Phe"] (D-loop context). Other metazoan
       groups may differ — consult the literature and set refinement.control_region_flanks
       accordingly. If a flank pair locates the gap, confidence is "high" regardless of AT%.

    2. Tandem repeat signal (secondary, discriminating):
       Among AT-rich candidate gaps (AT >= min_at, length >= min_len), gaps with detectable
       tandem repeats (~block_size bp units, >= min_copies copies at >= min_identity) are
       preferred. Confidence is "high" when a repeat is detected.

    3. AT% alone (fallback):
       When neither flanks nor repeats give a high-confidence result, the most AT-rich
       gap passing the length/AT filter is selected with confidence "low" (requires review).
    """
    gaps = _intergenic_gaps(record)
    if not gaps:
        return None

    # Stage 1: tRNA-flank detection (primary signal)
    flank_gap = None
    flank_method = None
    if control_region_flanks:
        flank_gap, flank_method = _find_control_region_by_flanks(record, control_region_flanks, gaps)

    if flank_gap is not None:
        best = flank_gap
        # Run tandem repeat as supporting evidence on the flank-identified gap.
        seq = (str(record.seq[best["start0"]:best["end0"]]) if not best["wrap"]
               else str(record.seq[best["start0"]:] + record.seq[:best["end0"]]))
        has_repeat, repeat_stats = _tandem_repeat_score(
            seq, repeat_block_size, repeat_min_identity, repeat_min_copies
        )
        best["has_repeat"] = has_repeat
        best["repeat_max_copies"] = repeat_stats.get("max_copies_found", 0)
        confidence = "high"
        repeat_note = " + tandem repeat signal" if has_repeat else " (no tandem repeat in assembly — likely collapsed)"
        inference = f"MitoCurator: {flank_method}{repeat_note}"
        note = "putative mitochondrial control region"
    else:
        # Stage 2+3: AT filter + tandem repeat fallback
        candidates = [g for g in gaps if g["length"] >= min_len and g["at_content"] >= min_at]
        if not candidates:
            return None

        for gap in candidates:
            seq = (str(record.seq[gap["start0"]:gap["end0"]]) if not gap["wrap"]
                   else str(record.seq[gap["start0"]:] + record.seq[:gap["end0"]]))
            has_repeat, repeat_stats = _tandem_repeat_score(
                seq, repeat_block_size, repeat_min_identity, repeat_min_copies
            )
            gap["has_repeat"] = has_repeat
            gap["repeat_max_copies"] = repeat_stats.get("max_copies_found", 0)

        high_conf = [g for g in candidates if g["has_repeat"]]
        pool = high_conf if high_conf else candidates
        best = max(pool, key=lambda g: g["at_content"])
        confidence = "high" if best["has_repeat"] else "low"

        if confidence == "high":
            note = "putative mitochondrial control region"
            inference = "MitoCurator: intergenic AT-rich scan + tandem repeat signal"
        else:
            note = ("putative AT-rich region; control region candidate — "
                    "no tandem repeat detected, requires manual review")
            inference = "MitoCurator: intergenic AT-rich scan only"

    s, e, n = best["start0"], best["end0"], len(record.seq)
    loc = (CompoundLocation([FeatureLocation(s, n, strand=1), FeatureLocation(0, e, strand=1)])
           if best["wrap"] else FeatureLocation(s, e, strand=1))

    record.features.append(SeqFeature(
        location=loc, type="misc_feature",
        qualifiers={"note": [note], "inference": [inference]},
    ))
    return {
        "seqid":             record.id,
        "start":             s + 1,
        "end":               e,
        "length":            best["length"],
        "at_content":        best["at_content"],
        "feature_type":      "misc_feature",
        "note":              note,
        "confidence":        confidence,
        "repeat_max_copies": best.get("repeat_max_copies", 0),
        "between":           "intergenic_scan",
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


def find_missing_cds_candidates(config, record, expected_gene_tsv, out_tsv):
    min_nt = int(safe_get(config, ["refinement", "orf_min_nt"], 150))
    genetic_code = int(safe_get(config, ["project", "genetic_code"], 5))
    exp = pd.read_csv(expected_gene_tsv, sep="\t")
    missing = exp[(exp["type"] == "CDS") & (exp["status"] == "MISSING")]["gene"].tolist()
    gaps = sorted(_intergenic_gaps(record), key=lambda g: g["length"], reverse=True)
    rows = []
    cid = 1

    for gene in missing:
        top_for_gene = []
        for gap in gaps:
            seq = (record.seq[gap["start0"]:gap["end0"]] if not gap["wrap"] else (record.seq[gap["start0"]:] + record.seq[:gap["end0"]]))
            for strand, nuc in [("+", str(seq).upper()), ("-", str(seq.reverse_complement()).upper())]:
                for frame in (0, 1, 2):
                    i = frame
                    while i + 3 <= len(nuc):
                        codon = nuc[i:i+3]
                        if codon in STOP_CODONS_TABLE5:
                            i += 3
                            continue
                        j = i
                        stops = 0
                        while j + 3 <= len(nuc):
                            c2 = nuc[j:j+3]
                            if c2 in STOP_CODONS_TABLE5:
                                stops += 1
                                break
                            j += 3
                        orf_len = j - i
                        if orf_len >= min_nt:
                            aa_len, internal, term = _translate_stop_metrics(nuc[i:j], code=genetic_code)
                            top_for_gene.append((internal, -orf_len, strand, frame, i, j, aa_len, term, gap))
                        i = j + 3
        top_for_gene.sort(key=lambda x: (x[0], x[1]))
        for cand in top_for_gene[:5]:
            internal, neglen, strand, frame, i, j, aa_len, term, gap = cand
            l = -neglen
            if strand == "+":
                start0 = (gap["start0"] + i) % len(record.seq)
                end0 = (gap["start0"] + j) % len(record.seq)
            else:
                start0 = (gap["start0"] + i) % len(record.seq)
                end0 = (gap["start0"] + j) % len(record.seq)
            left, right = _context_features(record, min(start0, end0), max(start0, end0))
            context = f"gap_{gap['start0']+1}..{gap['end0']}|len={gap['length']}|AT={gap['at_content']}"
            rows.append({
                "gene": gene, "candidate_id": f"{gene}_cand{cid}", "seqid": record.id,
                "start": start0 + 1, "end": end0, "strand": strand, "length_nt": l, "length_aa": aa_len,
                "frame": frame, "internal_stop_count": internal, "terminal_stop": term,
                "overlaps_existing_feature": "no", "nearest_left_feature": left, "nearest_right_feature": right,
                "region_context": context,
                "decision_hint": "STRONG_CANDIDATE_REGION" if gap["length"] >= 500 else "LOW_PRIORITY",
                "comment": "ORF scan on intergenic regions only; no automatic rescue applied",
            })
            cid += 1

    cols = ["gene", "candidate_id", "seqid", "start", "end", "strand", "length_nt", "length_aa", "frame", "internal_stop_count", "terminal_stop", "overlaps_existing_feature", "nearest_left_feature", "nearest_right_feature", "region_context", "decision_hint", "comment"]
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="\t", index=False)


def _load_ref_proteins(ref_gb: Path, genetic_code: int) -> dict[str, str]:
    """Load reference protein sequences keyed by normalised gene name."""
    proteins: dict[str, str] = {}
    if not ref_gb.exists():
        return proteins
    for rec in SeqIO.parse(str(ref_gb), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            raw = (feat.qualifiers.get("gene") or feat.qualifiers.get("product") or [None])[0]
            if raw is None:
                continue
            gene = _normalize_token(raw)
            try:
                proteins[gene] = str(
                    Seq(str(feat.extract(rec.seq)).upper()).translate(table=genetic_code, to_stop=False)
                ).rstrip("*")
            except Exception:
                continue
    return proteins


def _prot_coverage(nt: str, ref_aa: str, genetic_code: int) -> float | None:
    """Global alignment coverage of ref_aa by the frame-0 translation of nt."""
    try:
        trimmed = nt[: (len(nt) // 3) * 3]
        query_aa = str(Seq(trimmed).translate(table=genetic_code, to_stop=False)).rstrip("*")
        if not query_aa or not ref_aa:
            return None
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        aln = next(aligner.align(ref_aa, query_aa))
        aligned_ref = sum(e - s for s, e in aln.aligned[0])
        return aligned_ref / len(ref_aa)
    except Exception:
        return None


def find_cds_refinement_candidates(config, record, problematic_features_tsv, out_tsv):
    window = int(safe_get(config, ["refinement", "cds_refinement_window"], 300))
    code = int(safe_get(config, ["project", "genetic_code"], 5))
    completeness_min_ratio = float(safe_get(config, ["local_consensus", "completeness_min_ratio"], 0.80))
    completeness_min_cov   = float(safe_get(config, ["local_consensus", "completeness_min_cov"],   0.70))

    ref_gb_path = safe_get(config, ["input", "reference_gb"], None)
    ref_proteins = _load_ref_proteins(Path(ref_gb_path), code) if ref_gb_path else {}

    rows = []
    cds_feats = [f for f in record.features if f.type == "CDS"]

    for feat in cds_feats:
        gene = _normalize_token(feat.qualifiers.get("gene", [feature_name(feat)])[0])
        start0, end0 = int(feat.location.start), int(feat.location.end)
        strand = feat.location.strand or 1
        nt = str(feat.extract(record.seq)).upper()
        old_aa_len, old_internal, _ = _translate_stop_metrics(nt, code=code)

        if old_internal > 0:
            # Existing logic: small window scan for a better boundary (±30 bp).
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
                rows.append({
                    "gene": gene, "old_start": start0 + 1, "old_end": end0,
                    "old_strand": "+" if strand == 1 else "-", "old_length_nt": end0 - start0,
                    "old_internal_stop_count": old_internal,
                    "candidate_start": ".", "candidate_end": ".", "candidate_strand": ".",
                    "candidate_length_nt": ".", "candidate_length_aa": ".", "candidate_frame": ".",
                    "candidate_internal_stop_count": ".", "candidate_terminal_stop": ".",
                    "delta_start": ".", "delta_end": ".",
                    "decision_hint": "NO_BETTER_CANDIDATE",
                    "comment": "No alternative candidate in window search",
                    "problem_reason": "INTERNAL_STOP",
                })
                continue

            _, ns, ne, aa_len, internal, terminal, ds, de = best
            hint = "NO_BETTER_CANDIDATE"
            if internal < old_internal:
                hint = "SUGGEST_REVIEW"
            if internal == 0 and abs((ne - ns) - (end0 - start0)) <= window:
                hint = "STRONG_CANDIDATE"
            rows.append({
                "gene": gene, "old_start": start0 + 1, "old_end": end0,
                "old_strand": "+" if strand == 1 else "-", "old_length_nt": end0 - start0,
                "old_internal_stop_count": old_internal,
                "candidate_start": ns + 1, "candidate_end": ne,
                "candidate_strand": "+" if strand == 1 else "-",
                "candidate_length_nt": ne - ns, "candidate_length_aa": aa_len,
                "candidate_frame": (ns % 3), "candidate_internal_stop_count": internal,
                "candidate_terminal_stop": terminal, "delta_start": ds, "delta_end": de,
                "decision_hint": hint,
                "comment": "Coordinate-only candidate scan; no GenBank changes applied",
                "problem_reason": "INTERNAL_STOP",
            })

        else:
            # 0 internal stops — check completeness against the reference.
            if gene not in ref_proteins:
                continue
            ref_aa = ref_proteins[gene]
            ref_nt_len = len(ref_aa) * 3
            if ref_nt_len == 0:
                continue
            length_ratio = len(nt) / ref_nt_len
            cov = _prot_coverage(nt, ref_aa, code)

            problem_reason = None
            if length_ratio < completeness_min_ratio:
                problem_reason = "INCOMPLETE_LENGTH"
            elif cov is not None and cov < completeness_min_cov:
                problem_reason = "INCOMPLETE_COVERAGE"

            if problem_reason is None:
                continue

            comment = f"completeness: length_ratio={length_ratio:.3f}"
            if cov is not None:
                comment += f" prot_cov={cov:.3f}"
            comment += "; tblastn step in local_consensus will refine coordinates"

            rows.append({
                "gene": gene, "old_start": start0 + 1, "old_end": end0,
                "old_strand": "+" if strand == 1 else "-", "old_length_nt": end0 - start0,
                "old_internal_stop_count": 0,
                "candidate_start": start0 + 1, "candidate_end": end0,
                "candidate_strand": "+" if strand == 1 else "-",
                "candidate_length_nt": end0 - start0, "candidate_length_aa": old_aa_len,
                "candidate_frame": 0, "candidate_internal_stop_count": 0,
                "candidate_terminal_stop": ".", "delta_start": 0, "delta_end": 0,
                "decision_hint": "INCOMPLETE",
                "comment": comment,
                "problem_reason": problem_reason,
            })

    cols = [
        "gene", "old_start", "old_end", "old_strand", "old_length_nt",
        "old_internal_stop_count", "candidate_start", "candidate_end", "candidate_strand",
        "candidate_length_nt", "candidate_length_aa", "candidate_frame",
        "candidate_internal_stop_count", "candidate_terminal_stop",
        "delta_start", "delta_end", "decision_hint", "comment", "problem_reason",
    ]
    pd.DataFrame(rows, columns=cols).to_csv(out_tsv, sep="\t", index=False)


def refine_annotation(config, input_gb, outdir):
    outdir = ensure_dir(outdir)
    record, _ = read_record(input_gb)

    expected_tsv = outdir / "expected_gene_set.tsv"
    summarize_expected_gene_set(record, expected_tsv)

    min_len = int(safe_get(config, ["refinement", "at_rich_min_len"], 500))
    min_at = float(safe_get(config, ["refinement", "at_rich_min_at"], 75.0))
    annotate_at = bool(safe_get(config, ["refinement", "annotate_at_rich"], True))
    repeat_block_size = int(safe_get(config, ["refinement", "at_rich_repeat_block_size"], 55))
    repeat_min_identity = float(safe_get(config, ["refinement", "at_rich_repeat_min_identity"], 0.75))
    repeat_min_copies = int(safe_get(config, ["refinement", "at_rich_repeat_min_copies"], 3))
    # Default: Hymenoptera pair. Vertebrates: [["tRNA-Pro", "tRNA-Phe"]]. Set to [] to disable.
    control_region_flanks = safe_get(config, ["refinement", "control_region_flanks"],
                                     [["tRNA-Ile", "tRNA-Gln"]])
    added = add_at_rich_region(
        record, min_len=min_len, min_at=min_at,
        control_region_flanks=control_region_flanks,
        repeat_block_size=repeat_block_size,
        repeat_min_identity=repeat_min_identity,
        repeat_min_copies=repeat_min_copies,
    ) if annotate_at else None

    with open(outdir / "added_features.tsv", "w", encoding="utf-8") as out:
        out.write("seqid\tstart\tend\tlength\tat_content\tfeature_type\tnote\tconfidence\trepeat_max_copies\tbetween\n")
        if added:
            out.write(
                "{seqid}\t{start}\t{end}\t{length}\t{at_content}\t{feature_type}\t"
                "{note}\t{confidence}\t{repeat_max_copies}\t{between}\n".format(**added)
            )

    if bool(safe_get(config, ["refinement", "find_missing_cds_candidates"], True)):
        find_missing_cds_candidates(config, record, expected_tsv, outdir / "missing_gene_candidates.tsv")
    else:
        pd.DataFrame(columns=["gene", "candidate_id", "seqid", "start", "end", "strand", "length_nt", "length_aa", "frame", "internal_stop_count", "terminal_stop", "overlaps_existing_feature", "nearest_left_feature", "nearest_right_feature", "region_context", "decision_hint", "comment"]).to_csv(outdir / "missing_gene_candidates.tsv", sep="\t", index=False)

    find_cds_refinement_candidates(config, record, outdir / "problematic_features.tsv", outdir / "cds_refinement_candidates.tsv") if bool(safe_get(config, ["refinement", "find_cds_refinement_candidates"], True)) else pd.DataFrame(columns=["gene", "old_start", "old_end", "old_strand", "old_length_nt", "old_internal_stop_count", "candidate_start", "candidate_end", "candidate_strand", "candidate_length_nt", "candidate_length_aa", "candidate_frame", "candidate_internal_stop_count", "candidate_terminal_stop", "delta_start", "delta_end", "decision_hint", "comment"]).to_csv(outdir / "cds_refinement_candidates.tsv", sep="\t", index=False)

    out_gb = outdir / "refined.gb"
    write_record(record, out_gb, "genbank")
    return out_gb
