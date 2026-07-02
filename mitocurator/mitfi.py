from __future__ import annotations
from pathlib import Path
import datetime
import re
import subprocess

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

_AA_TO_THREE = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys',
    'L': 'Leu', 'L1': 'Leu', 'L2': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'S1': 'Ser', 'S2': 'Ser',
    'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
}


def _parse_mitfi_fasta(output: str) -> list[dict]:
    hits = []
    for line in output.splitlines():
        if not line.startswith('>'):
            continue
        parts = line[1:].split('|')
        if len(parts) < 9:
            continue
        hits.append({
            'seqid':     parts[0],
            'start':     int(parts[1]),    # 1-based
            'end':       int(parts[2]),    # 1-based inclusive
            'score':     float(parts[3]),
            'evalue':    parts[4],
            'anticodon': parts[5],
            'aa':        parts[6],
            'cm':        parts[7],
            'strand':    parts[8].strip(),
        })
    return hits


def _hit_to_feature(hit: dict) -> SeqFeature:
    three = _AA_TO_THREE.get(hit['aa'], hit['aa'])
    product = f"tRNA-{three}"
    strand = 1 if hit['strand'] == 'plus' else -1
    # BioPython FeatureLocation uses 0-based half-open [start, end)
    loc = FeatureLocation(hit['start'] - 1, hit['end'], strand=strand)
    feat = SeqFeature(loc, type='tRNA')
    feat.qualifiers['gene'] = [product]
    feat.qualifiers['product'] = [product]
    feat.qualifiers['note'] = [
        f"MiTFi v0.1; anticodon={hit['anticodon']}; "
        f"score={hit['score']}; evalue={hit['evalue']}; cm={hit['cm']}"
    ]
    return feat


def _trna_is_stop(feat) -> bool:
    """True if this tRNA feature is a tRNA-Stop (ARWEN false positive, anticodon TCA)."""
    product = feat.qualifiers.get('product', [''])[0].lower()
    return 'stop' in product


def _trna_is_artifact(feat) -> bool:
    """True if this tRNA would be flagged TRNA_LENGTH_ARTIFACT by gene_qc (< 50 or > 200 nt)."""
    length = int(feat.location.end) - int(feat.location.start)
    return length < 50 or length > 200


def _feature_score(feat) -> float | None:
    """Numeric score from a tRNA feature's note qualifier (MiTFi hits only)."""
    note = feat.qualifiers.get('note', [''])[0]
    for token in note.split(';'):
        token = token.strip()
        if token.startswith('score='):
            try:
                return float(token.split('=', 1)[1])
            except ValueError:
                return None
    return None


_ARWEN_LINE_RE = re.compile(
    r'mtRNA-\S+\s+(c?)\[(\d+),\s*(\d+)\]\s+\d+\s+\(([A-Za-z]+)\)'
)


def _find_arwen_file(source_gb) -> Path | None:
    """Locate MitoFinder's raw ARWEN output sibling to *source_gb*.

    MitoFinder writes  <run_dir>/<run_name>_arwen/<contig>.arwen
    alongside          <run_dir>/<run_name>_MitoFinder_arwen_Final_Results/<contig>.gb
    The GenBank conversion MitoFinder performs on ARWEN's hits drops the
    anticodon column, so it has to be recovered from this raw file. Returns
    None if it can't be located (e.g. MitoFinder was run with a different
    tRNA finder) — callers must then skip anticodon-aware logic.
    """
    if source_gb is None:
        return None
    source_gb = Path(source_gb)
    run_dir = source_gb.parent.parent
    candidate = run_dir / f"{run_dir.name}_arwen" / f"{source_gb.stem}.arwen"
    return candidate if candidate.exists() else None


def _parse_arwen_file(path: Path) -> dict[tuple[int, int, int], str]:
    """Map (start0, end0, strand) -> lowercase anticodon from ARWEN's raw output.

    Origin-wrapping hits (1-based start > end) are skipped: they can't be
    matched against a linear GenBank feature location, and in practice they
    are already filtered out upstream as TRNA_LENGTH_ARTIFACT.
    """
    lookup: dict[tuple[int, int, int], str] = {}
    for line in path.read_text().splitlines():
        m = _ARWEN_LINE_RE.search(line)
        if not m:
            continue
        minus, start1, end1, anticodon = m.groups()
        start1, end1 = int(start1), int(end1)
        if start1 > end1:
            continue
        strand = -1 if minus == 'c' else 1
        lookup[(start1 - 1, end1, strand)] = anticodon.lower()
    return lookup


def _attach_arwen_anticodons(record, source_gb) -> int:
    """Annotate ARWEN-derived tRNA features in *record* with qualifiers['anticodon'].

    Matches raw ARWEN hits (sibling to *source_gb*) to tRNA features by exact
    (start, end, strand). Returns the number of features annotated; a feature
    left unmatched (or if the .arwen file can't be found at all) simply has
    no 'anticodon' qualifier — callers must treat that as "unknown", not
    "no duplicate".
    """
    arwen_file = _find_arwen_file(source_gb)
    if arwen_file is None:
        return 0
    lookup = _parse_arwen_file(arwen_file)
    n_matched = 0
    for feat in record.features:
        if feat.type != 'tRNA':
            continue
        key = (int(feat.location.start), int(feat.location.end), feat.location.strand)
        anticodon = lookup.get(key)
        if anticodon:
            feat.qualifiers['anticodon'] = [anticodon]
            n_matched += 1
    return n_matched


def _product_to_aa(product: str) -> str:
    """Extract three-letter amino acid from tRNA product name.

    'tRNA-Phe' → 'Phe', 'tRNA-Leu2' → 'Leu', 'tRNA-Stop' → 'Stop'.
    Returns empty string if product is not a recognisable tRNA name.
    """
    p = product.strip()
    if not p.lower().startswith('trna-'):
        return ''
    return p[5:].rstrip('0123456789')


def apply_mitfi_fallback(
    record,
    mitfi_jar: Path,
    genetic_code: int,
    outdir: Path,
    expected_trna_count: int = 22,
    source_gb=None,
) -> dict:
    """Merge ARWEN and MiTFi tRNA predictions.

    Logic:
    1. Build arwen_valid: keep ARWEN tRNAs that are neither tRNA-Stop nor
       TRNA_LENGTH_ARTIFACT (< 50 or > 200 nt). Remove the rest from the record.
    1b. Deduplicate ARWEN over-predictions: multiple tRNA calls for the same
       amino acid AND same anticodon (recovered from MitoFinder's raw ARWEN
       file, sibling to source_gb) are true duplicates — keep the best-scoring
       one, or the first by position. Amino acid alone is NOT enough to group
       on: Leu and Ser legitimately have two isoacceptor tRNAs each (distinct
       anticodons), so those must never be collapsed. A tRNA whose anticodon
       can't be recovered is left untouched by this step.
    2. Identify amino acids absent from arwen_valid.
    3. Run MiTFi; for each MiTFi hit whose amino acid is absent from arwen_valid,
       add the MiTFi feature.  ARWEN always takes priority over MiTFi for covered
       amino acids.
    4. Trigger is unchanged: tRNA-Stop present OR total count < expected_trna_count.

    Returns audit dict with keys: mitfi_triggered, trigger_reason, trnas_before,
    trnas_after, arwen_kept, arwen_removed (list), arwen_deduplicated (dict),
    dedup_skipped_no_anticodon (list), mitfi_added (list).
    """
    _attach_arwen_anticodons(record, source_gb)

    all_trnas = [f for f in record.features if f.type == 'tRNA']
    trnas_before = len(all_trnas)
    has_stop = any(_trna_is_stop(f) for f in all_trnas)
    below_expected = trnas_before < expected_trna_count

    triggered = has_stop or below_expected
    reasons = []
    if has_stop:
        reasons.append('tRNA-Stop (anticodon TCA) detected — ARWEN false positive')
    if below_expected:
        reasons.append(f'tRNA count {trnas_before} < expected {expected_trna_count}')

    audit: dict = {
        'timestamp':       datetime.datetime.utcnow().isoformat(),
        'step':            'mitfi_trna_fallback',
        'mitfi_triggered': triggered,
        'trigger_reason':  '; '.join(reasons) if reasons else None,
        'trnas_before':    trnas_before,
        'trnas_after':     trnas_before,
        'arwen_kept':      trnas_before,
        'arwen_removed':   [],
        'arwen_deduplicated': {'n': 0, 'amino_acids': [], 'removed': []},
        'dedup_skipped_no_anticodon': [],
        'mitfi_added':     [],
    }

    if not triggered:
        return audit

    # ── Step 1: identify ARWEN tRNAs to remove ───────────────────────────────
    arwen_removed_info: list[dict] = []
    products_to_remove: set[str] = set()

    for feat in all_trnas:
        product = feat.qualifiers.get('product', [''])[0]
        if _trna_is_stop(feat):
            products_to_remove.add(product)
            arwen_removed_info.append({'product': product, 'reason': 'tRNA-Stop'})
        elif _trna_is_artifact(feat):
            length = int(feat.location.end) - int(feat.location.start)
            products_to_remove.add(product)
            arwen_removed_info.append({
                'product': product, 'length_nt': length, 'reason': 'TRNA_LENGTH_ARTIFACT',
            })

    def _flagged(feat) -> bool:
        """True for tRNA features and their gene wrappers that should be removed."""
        if feat.type == 'tRNA':
            return feat.qualifiers.get('product', [''])[0] in products_to_remove
        if feat.type == 'gene':
            return (feat.qualifiers.get('gene', [''])[0] in products_to_remove
                    or feat.qualifiers.get('product', [''])[0] in products_to_remove)
        return False

    record.features = [f for f in record.features if not _flagged(f)]

    # Amino acids still covered after filtering (arwen_valid)
    arwen_valid_trnas = [f for f in record.features if f.type == 'tRNA']

    # ── Step 1b: deduplicate ARWEN over-predictions — multiple tRNA calls for
    # the SAME amino acid AND SAME anticodon (true duplicates, e.g. tRNA-Leu2
    # x4 all reading anticodon 'taa'). Grouping by amino acid alone would
    # wrongly collapse legitimate isoacceptor pairs (Leu-UUR/taa vs
    # Leu-CUN/tag, same for Ser) into one. A feature whose anticodon
    # couldn't be recovered is never removed here — left visibly duplicate
    # is safer than discarding the wrong copy. ─────────────────────────────
    by_aa_anticodon: dict[tuple[str, str], list] = {}
    skipped_no_anticodon: list[dict] = []
    for feat in arwen_valid_trnas:
        product = feat.qualifiers.get('product', [''])[0]
        aa = _product_to_aa(product)
        if not aa:
            continue
        anticodon = feat.qualifiers.get('anticodon', [None])[0]
        if anticodon is None:
            skipped_no_anticodon.append({
                'product': product,
                'start':   int(feat.location.start),
                'end':     int(feat.location.end),
            })
            continue
        by_aa_anticodon.setdefault((aa, anticodon), []).append(feat)

    dedup_removed_info: list[dict] = []
    dedup_spans: set[tuple[int, int, str]] = set()
    for (aa, anticodon), feats in by_aa_anticodon.items():
        if len(feats) <= 1:
            continue
        ordered = sorted(
            feats,
            key=lambda f: (
                _feature_score(f) is None,
                -(_feature_score(f) or 0.0),
                int(f.location.start),
            ),
        )
        for feat in ordered[1:]:
            product = feat.qualifiers.get('product', [''])[0]
            dedup_spans.add((int(feat.location.start), int(feat.location.end), product))
            dedup_removed_info.append({
                'product':   product,
                'anticodon': anticodon,
                'start':     int(feat.location.start),
                'end':       int(feat.location.end),
            })

    if dedup_spans:
        def _is_dedup_span(f) -> bool:
            if f.type not in ('tRNA', 'gene'):
                return False
            product = f.qualifiers.get('gene', f.qualifiers.get('product', ['']))[0]
            return (int(f.location.start), int(f.location.end), product) in dedup_spans

        record.features = [f for f in record.features if not _is_dedup_span(f)]
        arwen_valid_trnas = [f for f in record.features if f.type == 'tRNA']

    arwen_aas: set[str] = set()
    for feat in arwen_valid_trnas:
        aa = _product_to_aa(feat.qualifiers.get('product', [''])[0])
        if aa:
            arwen_aas.add(aa)

    # ── Step 2-3: run MiTFi, add only amino acids absent from arwen_valid ────
    fasta_path = outdir / 'assembly_for_mitfi.fasta'
    SeqIO.write(record, str(fasta_path), 'fasta')

    mitfi_dir = mitfi_jar.parent
    cmd = [
        'java', '-jar', str(mitfi_jar),
        '-fasta', '-top', '-bstrands',
        '-code', str(genetic_code),
        str(fasta_path),
    ]
    result = subprocess.run(
        cmd, cwd=str(mitfi_dir), capture_output=True, text=True, check=True
    )

    with open(outdir / 'mitfi_raw_output.txt', 'w', encoding='utf-8') as fh:
        fh.write(result.stdout)

    hits = _parse_mitfi_fasta(result.stdout)

    mitfi_added_info: list[dict] = []
    mitfi_seen_aas: set[str] = set()

    for hit in hits:
        aa = _AA_TO_THREE.get(hit['aa'], hit['aa'])
        if aa in arwen_aas or aa in mitfi_seen_aas:
            continue
        record.features.append(_hit_to_feature(hit))
        mitfi_seen_aas.add(aa)
        mitfi_added_info.append({
            'product':   f"tRNA-{aa}",
            'anticodon': hit['anticodon'],
            'evalue':    hit['evalue'],
            'position':  f"{hit['start']}..{hit['end']}",
        })

    record.features.sort(key=lambda f: int(f.location.start))

    trnas_after = sum(1 for f in record.features if f.type == 'tRNA')
    audit['trnas_after']         = trnas_after
    audit['arwen_kept']          = len(arwen_valid_trnas)
    audit['arwen_removed']       = arwen_removed_info
    audit['arwen_deduplicated']  = {
        'n':            len(dedup_removed_info),
        'amino_acids':  sorted({aa for (aa, ac), feats in by_aa_anticodon.items() if len(feats) > 1}),
        'removed':      dedup_removed_info,
    }
    audit['dedup_skipped_no_anticodon'] = skipped_no_anticodon
    audit['mitfi_added']         = mitfi_added_info

    return audit
