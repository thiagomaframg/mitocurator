from __future__ import annotations

import csv
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .utils import ensure_dir, safe_get


@dataclass
class SequenceEdit:
    edit_id: str
    gene: str
    evidence_source: str
    decision: str
    start_1based: int
    end_1based: int
    reference_seq: str
    replacement_seq: str
    edit_type: str
    priority: str
    comment: str
    applied: bool = False
    status: str = "pending"
    old_length: int = 0
    new_length: int = 0
    delta: int = 0


def _read_tsv(path: Path) -> List[Dict[str, str]]:
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


def _parse_major_codons(value: str) -> List[str]:
    """Parse strings like hifi:TGA;illumina:TGA."""
    codons: List[str] = []
    for item in str(value or "").split(";"):
        item = item.strip()
        if not item or ":" not in item:
            continue
        _, codon = item.split(":", 1)
        codon = codon.strip().upper()
        if codon and codon not in {".", "N", "NNN"}:
            codons.append(codon)
    return codons


def _consensus_replacement_codon(value: str) -> str | None:
    codons = _parse_major_codons(value)
    unique = sorted(set(codons))
    if len(unique) == 1:
        return unique[0]
    return None


def _edit_type(ref: str, repl: str) -> str:
    if len(ref) == len(repl):
        return "SUBSTITUTION"
    if len(ref) < len(repl):
        return "INSERTION_OR_REPLACEMENT"
    return "DELETION_OR_REPLACEMENT"


def _feature_gene_name(feature) -> str:
    q = feature.qualifiers
    for key in ["gene", "locus_tag", "product", "label"]:
        if key in q and q[key]:
            return str(q[key][0])
    return "."


def _gene_strand_map(record: SeqRecord) -> Dict[str, int]:
    strands: Dict[str, int] = {}
    for feature in record.features:
        if feature.type not in {"CDS", "gene"}:
            continue
        gene = _feature_gene_name(feature)
        if not gene or gene == ".":
            continue

        parts = getattr(feature.location, "parts", None)
        if parts:
            strand = parts[0].strand
        else:
            strand = feature.location.strand

        if strand in {-1, 1}:
            # Prefer CDS strand when both gene and CDS are present.
            if feature.type == "CDS" or gene not in strands:
                strands[gene] = int(strand)
    return strands


def _seq_for_genomic_strand(coding_seq: str, strand: int) -> str:
    coding_seq = str(coding_seq or "").upper()
    if strand == -1:
        return str(Seq(coding_seq).reverse_complement())
    return coding_seq


def _collect_sequence_edits(root: Path, record: SeqRecord | None = None) -> List[SequenceEdit]:
    """Collect high-confidence sequence edits from read-support evidence.

    Current automatic application policy:
    - apply only rows where all readsets support the same coding replacement
    - require consensus_recommendation == CORRECTION_SUPPORTED_BY_ALL_READSETS
    - require replacement codon to be unambiguous
    """
    read_support_tsv = root / "10_read_support" / "readset_consensus_recommendations.tsv"
    rows = _read_tsv(read_support_tsv)

    strand_by_gene = _gene_strand_map(record) if record is not None else {}

    edits: List[SequenceEdit] = []

    for row in rows:
        recommendation = str(row.get("consensus_recommendation", "") or "")
        if recommendation != "CORRECTION_SUPPORTED_BY_ALL_READSETS":
            continue

        repl = _consensus_replacement_codon(row.get("major_codon_by_read_set", ""))
        if not repl:
            continue

        gene = str(row.get("gene", ".") or ".")
        stop = str(row.get("stop_aa_position", ".") or ".")
        start = int(row["genomic_codon_start"])
        end = int(row["genomic_codon_end"])
        ref_coding = str(row.get("reference_codon", "") or "").upper()
        strand = int(strand_by_gene.get(gene, 1))

        # Coordinates in read-support are genomic coordinates, but codons in
        # readset_consensus_recommendations.tsv are reported in coding-strand
        # orientation. For reverse-strand CDSs, convert both reference and
        # replacement codons to genomic orientation before editing the sequence.
        ref_genomic = _seq_for_genomic_strand(ref_coding, strand)
        repl_genomic = _seq_for_genomic_strand(repl, strand)

        edits.append(
            SequenceEdit(
                edit_id=f"{gene}_stop{stop}_{start}_{end}",
                gene=gene,
                evidence_source=str(read_support_tsv),
                decision=recommendation,
                start_1based=start,
                end_1based=end,
                reference_seq=ref_genomic,
                replacement_seq=repl_genomic,
                edit_type=_edit_type(ref_genomic, repl_genomic),
                priority=str(row.get("priority", "HIGH") or "HIGH"),
                comment=(
                    f"Applied because all evaluated readsets support the same "
                    f"non-stop replacement codon. Coding-strand edit: "
                    f"{ref_coding}->{repl}; genomic-strand edit: "
                    f"{ref_genomic}->{repl_genomic}; strand={strand}."
                ),
            )
        )

    return edits


def _validate_and_sort_edits(seq: str, edits: List[SequenceEdit]) -> List[SequenceEdit]:
    seq_upper = seq.upper()
    sorted_edits = sorted(edits, key=lambda e: (e.start_1based, e.end_1based))
    valid: List[SequenceEdit] = []

    last_end = 0
    for edit in sorted_edits:
        edit.old_length = edit.end_1based - edit.start_1based + 1
        edit.new_length = len(edit.replacement_seq)
        edit.delta = edit.new_length - edit.old_length

        if edit.start_1based < 1 or edit.end_1based > len(seq):
            edit.status = "skipped_out_of_bounds"
            continue

        if edit.start_1based <= last_end:
            edit.status = "skipped_overlapping_edit"
            continue

        observed = seq_upper[edit.start_1based - 1 : edit.end_1based]
        expected = edit.reference_seq.upper()

        if expected and observed != expected:
            edit.status = f"skipped_reference_mismatch_observed_{observed}"
            continue

        edit.status = "validated"
        valid.append(edit)
        last_end = edit.end_1based

    return valid


def _apply_sequence_edits(seq: str, edits: List[SequenceEdit]) -> str:
    """Apply edits using original coordinates.

    Coordinates are 1-based inclusive and refer to the original sequence.
    """
    parts: List[str] = []
    cursor0 = 0

    for edit in edits:
        start0 = edit.start_1based - 1
        end0 = edit.end_1based

        parts.append(seq[cursor0:start0])
        parts.append(edit.replacement_seq)

        cursor0 = end0
        edit.applied = True
        edit.status = "applied"

    parts.append(seq[cursor0:])
    return "".join(parts)


def _lift_boundary(pos0: int, edits: List[SequenceEdit], is_end: bool = False) -> int:
    """Lift a 0-based feature boundary from original to edited sequence."""
    shift = 0

    for edit in edits:
        start0 = edit.start_1based - 1
        end0 = edit.end_1based
        delta = edit.delta

        if pos0 < start0:
            break

        if pos0 > end0 or pos0 == end0:
            shift += delta
            continue

        if start0 <= pos0 < end0:
            if is_end:
                return start0 + shift + len(edit.replacement_seq)
            return start0 + shift

    return pos0 + shift


def _lift_location(location, edits: List[SequenceEdit]):
    parts = getattr(location, "parts", None)

    if parts:
        lifted_parts = []
        for part in parts:
            new_start = _lift_boundary(int(part.start), edits, is_end=False)
            new_end = _lift_boundary(int(part.end), edits, is_end=True)
            if new_end < new_start:
                new_end = new_start
            lifted_parts.append(FeatureLocation(new_start, new_end, strand=part.strand))

        # Biopython requires CompoundLocation to have at least two parts.
        # Some records may encode single-part features as a location with
        # a .parts attribute, so after liftover we safely collapse them back
        # to a simple FeatureLocation.
        if len(lifted_parts) == 1:
            return lifted_parts[0]

        return CompoundLocation(lifted_parts, operator=getattr(location, "operator", "join"))

    new_start = _lift_boundary(int(location.start), edits, is_end=False)
    new_end = _lift_boundary(int(location.end), edits, is_end=True)
    if new_end < new_start:
        new_end = new_start
    return FeatureLocation(new_start, new_end, strand=location.strand)


def _feature_name(feature: SeqFeature) -> str:
    q = feature.qualifiers
    for key in ["gene", "product", "locus_tag", "label"]:
        if key in q and q[key]:
            return str(q[key][0])
    return feature.type


def _feature_bounds(location) -> Tuple[int, int]:
    parts = getattr(location, "parts", None)
    if parts:
        starts = [int(p.start) for p in parts]
        ends = [int(p.end) for p in parts]
        return min(starts) + 1, max(ends)
    return int(location.start) + 1, int(location.end)


def _update_cds_translation(
    feature: SeqFeature,
    record: SeqRecord,
    genetic_code: int,
) -> SeqFeature:
    if feature.type != "CDS":
        return feature

    qualifiers = deepcopy(feature.qualifiers)
    qualifiers.pop("translation", None)

    try:
        nt = feature.extract(record.seq)
        if len(nt) % 3 == 0:
            aa = str(nt.translate(table=genetic_code, to_stop=False))
            if aa.endswith("*"):
                aa = aa[:-1]
            if "*" not in aa:
                qualifiers["translation"] = [aa]
            else:
                notes = qualifiers.get("note", [])
                notes.append("MitoCurator: translation not updated because internal stop remains after applied curation.")
                qualifiers["note"] = notes
        else:
            notes = qualifiers.get("note", [])
            notes.append("MitoCurator: translation not updated because CDS length is not a multiple of 3 after applied curation.")
            qualifiers["note"] = notes
    except Exception as exc:
        notes = qualifiers.get("note", [])
        notes.append(f"MitoCurator: translation not updated after applied curation ({exc}).")
        qualifiers["note"] = notes

    return SeqFeature(
        location=feature.location,
        type=feature.type,
        id=feature.id,
        qualifiers=qualifiers,
    )


def _lift_features(
    original_record: SeqRecord,
    curated_record: SeqRecord,
    edits: List[SequenceEdit],
    genetic_code: int,
) -> Tuple[List[SeqFeature], List[Dict[str, Any]]]:
    lifted_features: List[SeqFeature] = []
    rows: List[Dict[str, Any]] = []

    for idx, feature in enumerate(original_record.features, start=1):
        old_start, old_end = _feature_bounds(feature.location)
        new_location = _lift_location(feature.location, edits)
        new_start, new_end = _feature_bounds(new_location)

        qualifiers = deepcopy(feature.qualifiers)
        notes = qualifiers.get("note", [])
        notes.append("MitoCurator: coordinates lifted after evidence-backed sequence curation.")
        qualifiers["note"] = notes

        lifted = SeqFeature(
            location=new_location,
            type=feature.type,
            id=feature.id,
            qualifiers=qualifiers,
        )

        lifted.location = new_location
        lifted = _update_cds_translation(lifted, curated_record, genetic_code)

        lifted_features.append(lifted)

        rows.append(
            {
                "feature_index": idx,
                "feature_type": feature.type,
                "feature_name": _feature_name(feature),
                "old_start": old_start,
                "old_end": old_end,
                "old_length": old_end - old_start + 1,
                "new_start": new_start,
                "new_end": new_end,
                "new_length": new_end - new_start + 1,
                "coordinate_delta_start": new_start - old_start,
                "coordinate_delta_end": new_end - old_end,
            }
        )

    return lifted_features, rows


def _validate_curated_cds(record: SeqRecord, genetic_code: int) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene = _feature_name(feature)
        start, end = _feature_bounds(feature.location)

        try:
            nt = feature.extract(record.seq)
            aa = nt.translate(table=genetic_code, to_stop=False)
            aa_str = str(aa)
            terminal_stop = "yes" if aa_str.endswith("*") else "no"
            internal_stop_count = aa_str[:-1].count("*") if aa_str else 0
            internal_stop_positions = [
                str(i + 1) for i, x in enumerate(aa_str[:-1]) if x == "*"
            ]
            length_multiple_of_3 = "yes" if len(nt) % 3 == 0 else "no"
            status = "OK"
            if len(nt) % 3 != 0:
                status = "CDS_LENGTH_NOT_MULTIPLE_OF_3"
            elif internal_stop_count > 0:
                status = "INTERNAL_STOP_REMAINS"
        except Exception as exc:
            nt = ""
            aa_str = ""
            terminal_stop = "."
            internal_stop_count = -1
            internal_stop_positions = ["."]
            length_multiple_of_3 = "."
            status = f"VALIDATION_ERROR:{exc}"

        rows.append({
            "gene": gene,
            "start": start,
            "end": end,
            "strand": getattr(feature.location, "strand", "."),
            "cds_length_nt": len(nt),
            "length_multiple_of_3": length_multiple_of_3,
            "aa_length": len(aa_str),
            "terminal_stop": terminal_stop,
            "internal_stop_count": internal_stop_count,
            "internal_stop_positions": ";".join(internal_stop_positions) if internal_stop_positions else ".",
            "post_curation_status": status,
        })

    return rows


def _write_report(
    path: Path,
    edits: List[SequenceEdit],
    liftover_rows: List[Dict[str, Any]],
    original_len: int,
    curated_len: int,
) -> None:
    applied = [e for e in edits if e.applied]
    skipped = [e for e in edits if not e.applied]

    lines = []
    lines.append("# MitoCurator applied curation report\n")
    lines.append("This report records evidence-backed sequence edits applied to the mitogenome and the resulting coordinate liftover of GenBank features.\n")
    lines.append("## Summary\n")
    lines.append(f"- Original sequence length: {original_len}")
    lines.append(f"- Curated sequence length: {curated_len}")
    lines.append(f"- Net length change: {curated_len - original_len}")
    lines.append(f"- Applied edits: {len(applied)}")
    lines.append(f"- Skipped edits: {len(skipped)}")
    lines.append(f"- Features lifted: {len(liftover_rows)}\n")

    lines.append("## Applied edits\n")
    if applied:
        for e in applied:
            lines.append(
                f"- `{e.edit_id}` ({e.gene}): {e.start_1based}-{e.end_1based} "
                f"{e.reference_seq} → {e.replacement_seq}; delta={e.delta}; evidence={e.decision}"
            )
    else:
        lines.append("- No edits applied.")
    lines.append("")

    if skipped:
        lines.append("## Skipped edits\n")
        for e in skipped:
            lines.append(f"- `{e.edit_id}` ({e.gene}): {e.status}")
        lines.append("")

    lines.append("## Coordinate liftover\n")
    lines.append("All feature coordinates were recalculated from the original GenBank coordinates using the cumulative edit deltas. Features downstream of insertions/deletions are shifted automatically.\n")

    path.write_text("\n".join(lines), encoding="utf-8")


def run_applied_curation(config: dict, root: Path, outdir: Path | None = None) -> Path:
    outdir = ensure_dir(outdir or root / "15_applied_curation")

    input_gb = root / "05_refinement" / "refined.gb"
    if not input_gb.exists():
        raise FileNotFoundError(f"Missing refined GenBank: {input_gb}")

    original_record = SeqIO.read(str(input_gb), "genbank")
    original_seq = str(original_record.seq)

    genetic_code = int(
        safe_get(
            config,
            ["applied_curation", "genetic_code"],
            safe_get(config, ["mitofinder", "organism_code"], 5),
        )
    )

    raw_edits = _collect_sequence_edits(root, original_record)
    valid_edits = _validate_and_sort_edits(original_seq, raw_edits)
    curated_seq = _apply_sequence_edits(original_seq, valid_edits)

    curated_record = deepcopy(original_record)
    curated_record.seq = Seq(curated_seq)
    curated_record.id = original_record.id
    curated_record.name = original_record.name
    curated_record.description = original_record.description + " | MitoCurator evidence-backed curated sequence"

    lifted_features, liftover_rows = _lift_features(
        original_record,
        curated_record,
        valid_edits,
        genetic_code,
    )
    curated_record.features = lifted_features

    # Preserve molecule_type for GenBank writing.
    curated_record.annotations = deepcopy(original_record.annotations)
    curated_record.annotations.setdefault("molecule_type", "DNA")

    edit_rows = []
    for e in raw_edits:
        edit_rows.append(
            {
                "edit_id": e.edit_id,
                "gene": e.gene,
                "evidence_source": e.evidence_source,
                "decision": e.decision,
                "start_1based": e.start_1based,
                "end_1based": e.end_1based,
                "reference_seq": e.reference_seq,
                "replacement_seq": e.replacement_seq,
                "edit_type": e.edit_type,
                "old_length": e.old_length,
                "new_length": e.new_length,
                "delta": e.delta,
                "priority": e.priority,
                "applied": "yes" if e.applied else "no",
                "status": e.status,
                "comment": e.comment,
            }
        )

    _write_tsv(
        outdir / "applied_edits.tsv",
        edit_rows,
        [
            "edit_id",
            "gene",
            "evidence_source",
            "decision",
            "start_1based",
            "end_1based",
            "reference_seq",
            "replacement_seq",
            "edit_type",
            "old_length",
            "new_length",
            "delta",
            "priority",
            "applied",
            "status",
            "comment",
        ],
    )

    _write_tsv(
        outdir / "feature_coordinate_liftover.tsv",
        liftover_rows,
        [
            "feature_index",
            "feature_type",
            "feature_name",
            "old_start",
            "old_end",
            "old_length",
            "new_start",
            "new_end",
            "new_length",
            "coordinate_delta_start",
            "coordinate_delta_end",
        ],
    )

    post_cds_rows = _validate_curated_cds(curated_record, genetic_code)
    _write_tsv(
        outdir / "post_curation_cds_validation.tsv",
        post_cds_rows,
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
            "post_curation_status",
        ],
    )

    SeqIO.write(curated_record, str(outdir / "curated_mitogenome.gb"), "genbank")
    SeqIO.write(curated_record, str(outdir / "curated_mitogenome.fasta"), "fasta")

    _write_report(
        outdir / "applied_curation_report.md",
        raw_edits,
        liftover_rows,
        len(original_seq),
        len(curated_seq),
    )

    return outdir
