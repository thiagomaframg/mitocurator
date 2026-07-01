from __future__ import annotations
from pathlib import Path
import datetime
import json
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


def _count_trnas(record) -> int:
    return sum(1 for f in record.features if f.type == 'tRNA')


def _has_trna_stop(record) -> bool:
    for feat in record.features:
        if feat.type != 'tRNA':
            continue
        product = feat.qualifiers.get('product', [''])[0].lower()
        if 'stop' in product:
            return True
    return False


def _remove_trna_features(record) -> None:
    def _is_trna(feat) -> bool:
        if feat.type == 'tRNA':
            return True
        if feat.type == 'gene':
            gene = feat.qualifiers.get('gene', [''])[0].lower()
            product = feat.qualifiers.get('product', [''])[0].lower()
            return gene.startswith('trna') or product.startswith('trna')
        return False
    record.features = [f for f in record.features if not _is_trna(f)]


def apply_mitfi_fallback(
    record,
    mitfi_jar: Path,
    genetic_code: int,
    outdir: Path,
    expected_trna_count: int = 22,
) -> dict:
    """Check trigger conditions; if met, replace tRNA features with MiTFi predictions.

    Returns audit dict with keys: mitfi_triggered, trigger_reason, trnas_before,
    trnas_after, trnas_replaced, timestamp.
    """
    trnas_before = _count_trnas(record)
    has_stop = _has_trna_stop(record)
    below_expected = trnas_before < expected_trna_count

    triggered = has_stop or below_expected
    reasons = []
    if has_stop:
        reasons.append('tRNA-Stop (anticodon TCA) detected — ARWEN false positive')
    if below_expected:
        reasons.append(f'tRNA count {trnas_before} < expected {expected_trna_count}')

    audit: dict = {
        'timestamp':      datetime.datetime.utcnow().isoformat(),
        'step':           'mitfi_trna_fallback',
        'mitfi_triggered': triggered,
        'trigger_reason': '; '.join(reasons) if reasons else None,
        'trnas_before':   trnas_before,
        'trnas_after':    trnas_before,
        'trnas_replaced': 0,
    }

    if not triggered:
        return audit

    # Write FASTA extracted from the annotated record for MiTFi input.
    # Using the record sequence guarantees coordinate compatibility with the GB annotations.
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
    new_features = [_hit_to_feature(h) for h in hits]

    _remove_trna_features(record)
    record.features.extend(new_features)
    record.features.sort(key=lambda f: int(f.location.start))

    audit['trnas_after'] = len(new_features)
    audit['trnas_replaced'] = trnas_before

    return audit
