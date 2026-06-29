"""Local consensus repair for CDS with internal stops or missing from annotation.

Flow per CDS (see docs/mitocurator_dev_brief.md):
  [0] tblastn of reference protein vs assembly → correct CDS coordinates
  [1] extract assembly region ± flank_bp → region.fa
  [2] minimap2 | samtools sort/filter/index → mapped.filtered.bam
  [3] extract mapped reads → reads_candidate.fa
  [4] tblastn reference protein vs reads_candidate.fa → hits_vs_reads.tsv
  [5] extract aligned region from each read with hit (coding orientation)
  [6] filter extracted regions: translate 3 frames, keep ≤ max_stops_in_read
  [7] MAFFT filtered regions → majority-vote consensus
  [8] boundary shift scan (±shift_nt, 3 frames) → final candidate
  [9] validate; audit log JSONL + candidate FASTA (suggest) or in-place edit (apply)
"""
from __future__ import annotations

import json
import subprocess
from collections import Counter
from datetime import datetime
from pathlib import Path

import pysam
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

from .refinement import _normalize_token
from .utils import safe_get

_MINIMAP2_PRESETS: dict[str, str] = {
    "hifi":     "-ax map-hifi",
    "ont":      "-ax map-ont",
    "clr":      "-ax map-pb",
    "illumina": "-ax sr",
}


# ── Public ─────────────────────────────────────────────────────────────────────

def repair_cds_local_consensus(
    config: dict,
    record,
    problems: list[dict],
    outdir: Path,
    mode: str = "suggest",
    audit_log: Path | None = None,
) -> list[dict]:
    """Repair each CDS in *problems* using read-based local consensus.

    Parameters
    ----------
    config:     loaded config.yaml dict
    record:     SeqRecord; modified in-place when mode="apply"
    problems:   list of problem dicts from find_cds_refinement_candidates() or
                find_missing_cds_candidates(); must have fields gene, start (1-based),
                end (0-based exclusive), strand, internal_stop_count (or old_* variants)
    outdir:     pipeline output directory
    mode:       "diagnose" | "suggest" (default) | "apply"
    audit_log:  path to JSONL; defaults to <outdir>/06_local_consensus/audit_log.jsonl
    Returns list of audit entries written this call.
    """
    outdir = Path(outdir)
    step_name = safe_get(config, ["output", "step_dirs", "local_consensus"], "06_local_consensus")
    step_dir = outdir / step_name
    step_dir.mkdir(parents=True, exist_ok=True)

    if audit_log is None:
        audit_log = step_dir / "audit_log.jsonl"

    lc = safe_get(config, ["local_consensus"], {})
    flank_bp          = int(  safe_get(lc,     ["flank_bp"],               200))
    min_mq            = int(  safe_get(lc,     ["min_mq"],                  20))
    maj_freq          = float(safe_get(lc,     ["maj_freq"],               0.70))
    gap_col_thr       = float(safe_get(lc,     ["mafft_gap_col_threshold"], 0.50))
    max_n_frac        = float(safe_get(lc,     ["max_n_fraction"],          0.05))
    threads           = int(  safe_get(lc,     ["threads"],                  8))
    tblastn_evalue    = float(safe_get(lc,     ["tblastn_evalue"],          1e-5))
    shift_nt          = int(  safe_get(lc,     ["shift_nt"],                60))
    length_tol_frac   = float(safe_get(lc,     ["length_tol_frac"],         0.15))
    max_stops_in_read = int(  safe_get(lc,     ["max_stops_in_read"],        0))
    min_identity_pct  = float(safe_get(lc,     ["min_identity_pct"],        70.0))
    genetic_code      = int(  safe_get(config, ["project", "genetic_code"],  5))

    reads_cfg   = safe_get(config, ["reads"], {})
    technology  = safe_get(reads_cfg, ["technology"], "hifi")
    preset      = _MINIMAP2_PRESETS.get(technology, "-x map-hifi")
    reads_paths = _resolve_reads(reads_cfg)

    ref_gb       = safe_get(config, ["input", "reference_gb"], None)
    ref_proteins = _load_ref_proteins(Path(ref_gb), genetic_code) if ref_gb else {}

    versions       = _tool_versions()
    audit_entries: list[dict] = []
    summary_rows:  list[dict] = []

    genome_fa = step_dir / "genome.fa"
    genome_fa.write_text(f">{record.id}\n{str(record.seq).upper()}\n")

    for p in problems:
        gene     = p["gene"]
        gene_dir = step_dir / gene
        gene_dir.mkdir(exist_ok=True)

        start1     = int(p.get("start", p.get("old_start", p.get("candidate_start", 1))))
        end0       = int(p.get("end",   p.get("old_end",   p.get("candidate_end", len(record.seq)))))
        start0     = start1 - 1
        strand_raw = p.get("strand", p.get("old_strand", p.get("candidate_strand", "+")))
        strand     = 1 if strand_raw in (1, "+") else -1
        stops_bef  = int(p.get("internal_stop_count", p.get("old_internal_stop_count", 0)))
        prob_type  = "PROBLEM_INTERNAL_STOP" if stops_bef > 0 else "MISSING"

        if gene not in ref_proteins:
            entry = _make_entry(gene, prob_type, "SKIPPED_NO_REFERENCE", versions, [],
                                stops_bef, stops_bef)
            _append_audit(audit_log, entry)
            audit_entries.append(entry)
            continue

        if not reads_paths:
            entry = _make_entry(gene, prob_type, "SKIPPED_NO_READS", versions, [],
                                stops_bef, stops_bef)
            _append_audit(audit_log, entry)
            audit_entries.append(entry)
            continue

        # ── [0] tblastn coordinate refinement ────────────────────────────────
        prot_fa     = gene_dir / f"{gene}_ref.fa"
        tblastn_tsv = gene_dir / "tblastn.tsv"
        prot_seq    = ref_proteins[gene].rstrip("*")
        prot_fa.write_text(f">{gene}\n{prot_seq}\n")
        tblastn_cmd = _run_tblastn(prot_fa, genome_fa, genetic_code, tblastn_tsv, threads, tblastn_evalue)
        coords = _best_tblastn_coords(tblastn_tsv, len(prot_seq))
        coords_confirmed = coords is not None
        if coords_confirmed:
            new_s0, new_e0, new_strand = coords
            start0 = max(0, new_s0)
            end0   = min(len(record.seq), new_e0)
            strand = new_strand

        # ── [1] extract assembly region ───────────────────────────────────────
        target_fa = gene_dir / "region.fa"
        reg_start0, reg_end0, cds_in_start, cds_in_end = _extract_region(
            record, start0, end0, flank_bp, target_fa, gene
        )
        bam_path     = gene_dir / "mapped.filtered.bam"
        ref_name     = f"region_{gene}"
        cur_coord    = f"{start0} {end0} {strand}"
        region_label = f"{record.id}:{reg_start0}-{reg_end0}"

        # ── [2] map reads to assembly region ─────────────────────────────────
        coord_cache = gene_dir / "coords.txt"
        if (bam_path.exists()
                and coord_cache.exists()
                and coord_cache.read_text().strip() == cur_coord):
            cmds_map = []
        else:
            if bam_path.exists():
                bam_path.unlink()
                Path(str(bam_path) + ".bai").unlink(missing_ok=True)
            cmds_map = _map_and_filter(reads_paths, preset, target_fa, bam_path, threads, min_mq)
            coord_cache.write_text(cur_coord)
        cmds = [tblastn_cmd] + cmds_map

        # ── [3] extract mapped reads → reads_candidate.fa ────────────────────
        reads_fa = gene_dir / "reads_candidate.fa"
        n_reads_fetched = _extract_reads_from_bam(bam_path, ref_name, reads_fa)

        if n_reads_fetched == 0:
            entry = _make_entry(gene, prob_type, "SKIPPED_NO_READS_IN_BAM", versions, cmds,
                                stops_bef, stops_bef,
                                evidence_extra={"region": region_label})
            _append_audit(audit_log, entry)
            audit_entries.append(entry)
            continue

        # ── [4] tblastn protein vs reads ─────────────────────────────────────
        hits_tsv = gene_dir / "hits_vs_reads.tsv"
        tblastn_reads_cmd = _run_tblastn_vs_reads(
            prot_fa, reads_fa, genetic_code, hits_tsv, tblastn_evalue
        )
        cmds.append(tblastn_reads_cmd)

        # ── [5] extract aligned regions from reads ────────────────────────────
        hit_regions = _extract_hit_regions(reads_fa, hits_tsv, len(prot_seq))
        n_hits = len(hit_regions)

        if n_hits == 0:
            entry = _make_entry(gene, prob_type, "SKIPPED_NO_TBLASTN_HIT_IN_READS", versions, cmds,
                                stops_bef, stops_bef,
                                evidence_extra={
                                    "region": region_label,
                                    "reads_fetched_from_bam": n_reads_fetched,
                                    "reads_with_tblastn_hit": 0,
                                })
            _append_audit(audit_log, entry)
            audit_entries.append(entry)
            continue

        # ── [6] filter by internal stops (3 frames) ───────────────────────────
        filtered_seqs = _filter_reads_by_stops(hit_regions, genetic_code, max_stops_in_read)
        n_filtered = len(filtered_seqs)

        if n_filtered == 0:
            entry = _make_entry(gene, prob_type, "SKIPPED_NO_READS_PASS_FILTER", versions, cmds,
                                stops_bef, stops_bef,
                                evidence_extra={
                                    "region": region_label,
                                    "reads_fetched_from_bam": n_reads_fetched,
                                    "reads_with_tblastn_hit": n_hits,
                                    "reads_passing_stop_filter": 0,
                                    "max_stops_in_read": max_stops_in_read,
                                })
            _append_audit(audit_log, entry)
            audit_entries.append(entry)
            continue

        # ── [7] MAFFT + majority consensus ────────────────────────────────────
        consensus_nt = _mafft_consensus(filtered_seqs, gene_dir, threads, gap_col_thr, maj_freq)

        if not consensus_nt:
            entry = _make_entry(gene, prob_type, "SKIPPED_NO_CONSENSUS", versions, cmds,
                                stops_bef, stops_bef,
                                evidence_extra={
                                    "region": region_label,
                                    "reads_fetched_from_bam": n_reads_fetched,
                                    "reads_with_tblastn_hit": n_hits,
                                    "reads_passing_stop_filter": n_filtered,
                                })
            _append_audit(audit_log, entry)
            audit_entries.append(entry)
            continue

        # ── [8] boundary shift scan ───────────────────────────────────────────
        ref_prot_len = len(prot_seq)
        best_seq = _scan_boundary_shifts(
            consensus_nt, ref_prot_len, genetic_code, shift_nt, length_tol_frac
        )

        if best_seq is None:
            val_raw = _validate_candidate(consensus_nt, genetic_code, ref_proteins[gene])
            entry = _make_entry(gene, prob_type, "REJECTED_NO_BOUNDARY_FIT", versions, cmds,
                                stops_bef, val_raw["internal_stops"],
                                evidence_extra={
                                    "region": region_label,
                                    "reads_fetched_from_bam": n_reads_fetched,
                                    "reads_with_tblastn_hit": n_hits,
                                    "reads_passing_stop_filter": n_filtered,
                                    "max_stops_in_read": max_stops_in_read,
                                    "mafft_gap_col_threshold": gap_col_thr,
                                    "consensus_n_fraction": val_raw["n_fraction"],
                                    "ref_protein_identity_pct": val_raw.get("identity_pct"),
                                    "ref_protein_coverage_pct": val_raw.get("coverage_pct"),
                                })
            _append_audit(audit_log, entry)
            audit_entries.append(entry)
            continue

        # ── [9] validate + write + apply ──────────────────────────────────────
        val       = _validate_candidate(best_seq, genetic_code, ref_proteins[gene])
        stops_aft = val["internal_stops"]
        n_frac    = val["n_fraction"]

        identity_pct = val.get("identity_pct")

        if stops_aft > 0:
            action = "REJECTED_STOPS_REMAIN"
        elif n_frac > max_n_frac:
            action = "REJECTED_HIGH_N"
        elif identity_pct is not None and identity_pct < min_identity_pct:
            action = "REJECTED_LOW_IDENTITY"
        elif not coords_confirmed:
            # Candidate sequence is valid but assembly coordinates were not confirmed
            # by tblastn — splice position is unreliable; do not modify the record.
            action = "REJECTED_TBLASTN_NO_HIT_IN_ASSEMBLY"
        elif mode == "apply":
            action = "APPLIED"
        else:
            action = "SUGGEST"

        # Write candidate FASTA for all non-diagnose outcomes where a valid sequence
        # was built (including REJECTED_TBLASTN_NO_HIT_IN_ASSEMBLY for inspection).
        _will_write = action in ("SUGGEST", "APPLIED",
                                  "REJECTED_LOW_IDENTITY", "REJECTED_TBLASTN_NO_HIT_IN_ASSEMBLY")
        cand_fa = None
        if mode != "diagnose" and _will_write:
            cand_fa = gene_dir / f"{gene}_candidate.fa"
            cand_fa.write_text(f">{gene}_local_consensus\n{best_seq}\n")
            prot_out = gene_dir / f"{gene}_candidate_protein.fa"
            prot_out.write_text(f">{gene}_local_consensus_protein\n{val['aa_seq']}\n")

        if action == "APPLIED":
            _apply_to_record(record, gene, start0, end0, strand, best_seq, prob_type, genetic_code)

        entry = _make_entry(
            gene, prob_type, action, versions, cmds,
            stops_bef, stops_aft,
            evidence_extra={
                "region": region_label,
                "reads_fetched_from_bam": n_reads_fetched,
                "reads_with_tblastn_hit": n_hits,
                "reads_passing_stop_filter": n_filtered,
                "max_stops_in_read": max_stops_in_read,
                "mafft_gap_col_threshold": gap_col_thr,
                "consensus_n_fraction": n_frac,
                "assembly_coords_confirmed": coords_confirmed,
                "ref_protein_identity_pct": val.get("identity_pct"),
                "ref_protein_coverage_pct": val.get("coverage_pct"),
            },
            candidate={
                "sequence_nt": best_seq if _will_write else None,
                "length_nt": len(best_seq),
                "length_aa": len(val["aa_seq"]) - 1,
                "internal_stops": stops_aft,
                "terminal_stop": "yes" if val["terminal_stop"] else "no",
                "fasta_path": str(cand_fa) if cand_fa else None,
            },
        )
        _append_audit(audit_log, entry)
        audit_entries.append(entry)

        summary_rows.append({
            "gene": gene, "problem": prob_type, "action": action,
            "stops_before": stops_bef, "stops_after": stops_aft,
            "n_fraction": round(n_frac, 4), "length_nt": len(best_seq),
        })

    if summary_rows:
        import pandas as pd
        pd.DataFrame(summary_rows).to_csv(step_dir / "summary.tsv", sep="\t", index=False)

    return audit_entries


# ── Private helpers ────────────────────────────────────────────────────────────

def _resolve_reads(reads_cfg: dict) -> list[str]:
    paths: list[str] = []
    for key in ("hifi", "ont", "clr", "illumina_r1", "illumina_r2"):
        v = reads_cfg.get(key)
        if v is None:
            continue
        paths.extend(str(x) for x in (v if isinstance(v, list) else [v]))
    return [p for p in paths if Path(p).exists()]


def _load_ref_proteins(ref_gb: Path, genetic_code: int) -> dict[str, str]:
    proteins: dict[str, str] = {}
    if not ref_gb.exists():
        return proteins
    for rec in SeqIO.parse(str(ref_gb), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            raw = feat.qualifiers.get("gene", feat.qualifiers.get("product", [None]))[0]
            if raw is None:
                continue
            gene = _normalize_token(raw)
            try:
                aa = str(Seq(str(feat.extract(rec.seq)).upper()).translate(
                    table=genetic_code, to_stop=False))
            except Exception:
                continue
            proteins[gene] = aa
    return proteins


def _extract_region(record, start0: int, end0: int, flank_bp: int,
                    out_fa: Path, gene: str) -> tuple[int, int, int, int]:
    """Extract [start0-flank, end0+flank] from record, write FASTA.

    Returns (reg_start0, reg_end0, cds_start_in_region, cds_end_in_region).
    """
    seq_len    = len(record.seq)
    reg_start0 = max(0, start0 - flank_bp)
    reg_end0   = min(seq_len, end0 + flank_bp)
    left_flank = start0 - reg_start0
    cds_len    = end0 - start0
    out_fa.write_text(f">region_{gene}\n{str(record.seq[reg_start0:reg_end0]).upper()}\n")
    return reg_start0, reg_end0, left_flank, left_flank + cds_len


def _map_and_filter(reads_paths: list[str], preset: str, region_fa: Path,
                    bam_path: Path, threads: int, min_mq: int) -> list[str]:
    raw_bam  = bam_path.with_suffix(".unsorted.bam")
    map_cmd  = (f"minimap2 -t {threads} {preset} {region_fa} {' '.join(reads_paths)}"
                f" | samtools sort -@ {threads} -o {raw_bam}")
    filt_cmd = f"samtools view -b -q {min_mq} -F 2308 {raw_bam} -o {bam_path}"
    idx_cmd  = f"samtools index {bam_path}"
    subprocess.run(map_cmd,  shell=True, check=True)
    subprocess.run(filt_cmd, shell=True, check=True)
    subprocess.run(idx_cmd,  shell=True, check=True)
    raw_bam.unlink(missing_ok=True)
    return [map_cmd, filt_cmd, idx_cmd]


def _extract_reads_from_bam(bam_path: Path, ref_name: str, out_fa: Path) -> int:
    """Dump all non-duplicate read sequences from bam_path to out_fa FASTA.

    Deduplicates by query name (keeps first occurrence). Returns read count.
    """
    seen: set[str] = set()
    count = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        with open(out_fa, "w") as fa:
            for read in bam.fetch(ref_name):
                if read.is_unmapped or read.query_sequence is None:
                    continue
                if read.query_name in seen:
                    continue
                seen.add(read.query_name)
                fa.write(f">{read.query_name}\n{read.query_sequence}\n")
                count += 1
    return count


def _run_tblastn_vs_reads(prot_fa: Path, reads_fa: Path, genetic_code: int,
                           out_tsv: Path, evalue: float) -> str:
    # -num_threads is silently ignored by tblastn when -subject is specified (known limitation).
    cmd = (f"tblastn -query {prot_fa} -subject {reads_fa}"
           f" -db_gencode {genetic_code} -evalue {evalue}"
           f" -outfmt '6 qaccver saccver pident length qstart qend sstart send evalue bitscore'"
           f" -out {out_tsv}")
    subprocess.run(cmd, shell=True, check=True)
    return cmd


def _extract_hit_regions(reads_fa: Path, hits_tsv: Path, ref_prot_len: int) -> dict[str, str]:
    """Extract full expected CDS from each read with a tblastn hit.

    Uses qstart/qend to extend the tblastn hit window to cover the full expected
    CDS length (same extension logic as _best_tblastn_coords): partial hits still
    recover the complete gene region from the read. For minus-strand hits
    (sstart > send), returns the reverse complement.
    Returns {read_id: coding_subseq}.
    """
    reads: dict[str, str] = {
        rec.id: str(rec.seq).upper()
        for rec in SeqIO.parse(str(reads_fa), "fasta")
    }
    # best hit per read: (sstart, send, qstart, qend, bitscore)
    best: dict[str, tuple[int, int, int, int, float]] = {}
    try:
        with open(hits_tsv) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 10:
                    continue
                read_id  = parts[1]
                qstart   = int(parts[4])
                qend     = int(parts[5])
                sstart   = int(parts[6])
                send     = int(parts[7])
                bitscore = float(parts[9])
                if read_id not in best or bitscore > best[read_id][4]:
                    best[read_id] = (sstart, send, qstart, qend, bitscore)
    except (FileNotFoundError, ValueError):
        return {}

    regions: dict[str, str] = {}
    for read_id, (sstart, send, qstart, qend, _) in best.items():
        seq = reads.get(read_id)
        if seq is None:
            continue
        n          = len(seq)
        ext_5prime = (qstart - 1) * 3
        ext_3prime = (ref_prot_len - qend) * 3
        if sstart <= send:
            # + strand: protein N-term at sstart, C-term at send
            s0 = max(0, sstart - 1 - ext_5prime)
            e0 = min(n, send + ext_3prime)
            subseq = seq[s0:e0]
        else:
            # - strand: protein N-term at sstart (high coord), C-term at send (low coord)
            s0 = max(0, send - 1 - ext_3prime)
            e0 = min(n, sstart + ext_5prime)
            subseq = str(Seq(seq[s0:e0]).reverse_complement())
        if subseq:
            regions[read_id] = subseq
    return regions


def _filter_reads_by_stops(
    regions: dict[str, str], genetic_code: int, max_internal_stops: int
) -> list[str]:
    """Keep read regions where at least one reading frame has ≤ max_internal_stops.

    For each passing read, retains the subsequence in the best frame (fewest stops;
    ties broken by longest coding length).
    """
    result: list[str] = []
    for seq in regions.values():
        best_seq   = None
        best_stops = max_internal_stops + 1
        best_len   = 0
        for offset in range(3):
            sub     = seq[offset:]
            trimmed = sub[: (len(sub) // 3) * 3]
            if not trimmed:
                continue
            try:
                aa = str(Seq(trimmed).translate(table=genetic_code, to_stop=False))
            except Exception:
                continue
            internal = aa[:-1].count("*") if len(aa) > 1 else 0
            if internal < best_stops or (internal == best_stops and len(trimmed) > best_len):
                best_stops = internal
                best_seq   = trimmed
                best_len   = len(trimmed)
        if best_seq is not None and best_stops <= max_internal_stops:
            result.append(best_seq)
    return result


def _mafft_consensus(
    seqs: list[str], gene_dir: Path, threads: int, gap_col_thr: float, maj_freq: float
) -> str:
    """Write seqs to FASTA, run MAFFT, return column-majority consensus.

    Columns where gap_fraction > gap_col_thr are skipped (read-specific insertions
    relative to the majority; skipping prevents spurious bases from entering the consensus).
    """
    if len(seqs) == 1:
        return seqs[0]

    in_fa  = gene_dir / "mafft_input.fa"
    out_fa = gene_dir / "mafft_aligned.fa"

    with open(in_fa, "w") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f">read_{i}\n{seq}\n")

    cmd = f"mafft --auto --thread {threads} --quiet {in_fa} > {out_fa}"
    subprocess.run(cmd, shell=True, check=True)

    aligned = [str(rec.seq).upper() for rec in SeqIO.parse(str(out_fa), "fasta")]
    if not aligned:
        return ""

    n_seqs  = len(aligned)
    aln_len = max(len(s) for s in aligned)
    consensus: list[str] = []

    for col in range(aln_len):
        bases    = [s[col] for s in aligned if col < len(s)]
        n_gaps   = bases.count('-')
        gap_frac = n_gaps / n_seqs

        if gap_frac > gap_col_thr:
            continue

        non_gap = [b for b in bases if b != '-']
        if not non_gap:
            continue
        top_base, top_count = Counter(non_gap).most_common(1)[0]
        consensus.append(top_base if top_count / len(non_gap) >= maj_freq else 'N')

    return ''.join(consensus)


def _scan_boundary_shifts(
    consensus_nt: str,
    ref_prot_len: int,
    genetic_code: int,
    shift_nt: int,
    tol_frac: float,
) -> str | None:
    """Find best reading frame / boundary by trimming up to shift_nt nt from each end.

    Tests all (d5, d3) pairs where d5, d3 ∈ [0, shift_nt], and all 3 reading frames.
    Accepts candidates with:
      - 0 internal stop codons
      - |aa_len - ref_prot_len| / ref_prot_len ≤ tol_frac
    Returns the candidate with fewest stops then smallest length deviation, or None.
    """
    best_stops:    int       = 999
    best_len_diff: int       = 999
    best_seq:      str | None = None

    n = len(consensus_nt)

    for d5 in range(min(shift_nt + 1, n)):
        for d3 in range(min(shift_nt + 1, n - d5)):
            sub = consensus_nt[d5 : n - d3] if d3 > 0 else consensus_nt[d5:]
            if not sub:
                continue
            for frame in range(3):
                fseq    = sub[frame:]
                trimmed = fseq[: (len(fseq) // 3) * 3]
                if not trimmed:
                    continue
                try:
                    aa = str(Seq(trimmed).translate(table=genetic_code, to_stop=False))
                except Exception:
                    continue
                internal = aa[:-1].count("*") if len(aa) > 1 else 0
                aa_len   = len(aa.rstrip("*"))
                len_diff = abs(aa_len - ref_prot_len)
                if len_diff / ref_prot_len > tol_frac:
                    continue
                if internal == 0 and (
                    best_seq is None
                    or internal < best_stops
                    or (internal == best_stops and len_diff < best_len_diff)
                ):
                    best_stops    = internal
                    best_len_diff = len_diff
                    best_seq      = trimmed

    return best_seq


def _validate_candidate(seq_nt: str, genetic_code: int, ref_protein: str) -> dict:
    try:
        aa_full = str(Seq(seq_nt).translate(table=genetic_code, to_stop=False))
    except Exception:
        return {"internal_stops": 999, "terminal_stop": False,
                "n_fraction": 1.0, "aa_seq": "", "identity_pct": None, "coverage_pct": None}

    internal   = aa_full[:-1].count("*") if len(aa_full) > 1 else 0
    terminal   = bool(aa_full) and aa_full[-1] == "*"
    n_fraction = seq_nt.upper().count('N') / len(seq_nt) if seq_nt else 1.0

    identity_pct, coverage_pct = _protein_identity(ref_protein, aa_full)
    return {
        "internal_stops": internal,
        "terminal_stop":  terminal,
        "n_fraction":     round(n_fraction, 4),
        "aa_seq":         aa_full,
        "identity_pct":   identity_pct,
        "coverage_pct":   coverage_pct,
    }


def _protein_identity(ref_prot: str, query_prot: str) -> tuple[float | None, float | None]:
    ref_clean = ref_prot.rstrip("*")
    qry_clean = query_prot.rstrip("*")
    if not ref_clean or not qry_clean:
        return None, None
    try:
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score   = -11
        aligner.extend_gap_score = -1
        aln = next(aligner.align(ref_clean, qry_clean))
        ref_ranges, qry_ranges = aln.aligned
        matches = sum(
            sum(1 for a, b in zip(ref_clean[rs:re], qry_clean[qs:qe]) if a == b)
            for (rs, re), (qs, qe) in zip(ref_ranges, qry_ranges)
        )
        identity_pct = round(matches / len(ref_clean) * 100, 1)
        coverage_pct = round(len(qry_clean) / len(ref_clean) * 100, 1)
        return identity_pct, coverage_pct
    except Exception:
        return None, None


def _apply_to_record(record, gene: str, start0: int, end0: int, strand: int,
                     new_seq: str, prob_type: str, genetic_code: int) -> None:
    """Replace record.seq[start0:end0] with new_seq and shift all feature coordinates.

    new_seq must be in coding orientation (5'→3' of gene). For minus-strand genes
    the assembly stores the reverse complement, so we RC before splicing.
    """
    asm_seq = (str(Seq(new_seq).reverse_complement()) if strand == -1 else new_seq)
    old_len = end0 - start0
    new_len = len(asm_seq)
    delta   = new_len - old_len

    record.seq = record.seq[:start0] + Seq(asm_seq) + record.seq[end0:]

    target_gene    = _normalize_token(gene)
    target_updated = False

    for feat in record.features:
        fs      = int(feat.location.start)
        fe      = int(feat.location.end)
        fstrand = feat.location.strand

        if fe <= start0:
            continue

        if fs >= end0:
            feat.location = FeatureLocation(fs + delta, fe + delta, strand=fstrand)
            continue

        is_target = (_normalize_token(
            feat.qualifiers.get("gene", feat.qualifiers.get("product", [gene]))[0]
        ) == target_gene and feat.type == "CDS")

        if is_target:
            feat.location = FeatureLocation(start0, start0 + new_len, strand=fstrand)
            target_updated = True
        else:
            new_fe = (fe + delta) if fe > end0 else fe
            feat.location = FeatureLocation(fs, new_fe, strand=fstrand)

    if prob_type == "MISSING" and not target_updated:
        record.features.append(SeqFeature(
            FeatureLocation(start0, start0 + new_len, strand=strand),
            type="CDS",
            qualifiers={
                "gene":         [gene],
                "codon_start":  [1],
                "transl_table": [str(genetic_code)],
                "note":         ["APPLIED_local_consensus; recovered by read-based majority consensus"],
            },
        ))


def _run_tblastn(protein_fa: Path, genome_fa: Path, genetic_code: int,
                 out_tsv: Path, threads: int, evalue: float) -> str:
    cmd = (f"tblastn -query {protein_fa} -subject {genome_fa}"
           f" -db_gencode {genetic_code} -evalue {evalue}"
           f" -outfmt 6 -num_threads {threads} -out {out_tsv}")
    subprocess.run(cmd, shell=True, check=True)
    return cmd


def _best_tblastn_coords(tblastn_tsv: Path, query_len: int) -> tuple[int, int, int] | None:
    """Best tblastn format-6 hit → (start0, end0, strand).

    Extends hit by (qstart-1)*3 at 5' and (query_len-qend)*3 at 3' to recover
    unaligned protein termini. Returns 0-based half-open coordinates.
    """
    best_score = -1.0
    best: tuple[int, int, int, int] | None = None
    try:
        with open(tblastn_tsv) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 12:
                    continue
                qstart = int(parts[6])
                qend   = int(parts[7])
                sstart = int(parts[8])
                send   = int(parts[9])
                score  = float(parts[11])
                if score > best_score:
                    best_score = score
                    best = (sstart, send, qstart, qend)
    except (FileNotFoundError, ValueError):
        return None
    if best is None:
        return None

    sstart, send, qstart, qend = best
    ext_5prime = (qstart - 1) * 3
    ext_3prime = (query_len - qend) * 3

    if sstart <= send:
        start0 = sstart - 1 - ext_5prime
        end0   = send + ext_3prime
    else:
        start0 = send - 1 - ext_3prime
        end0   = sstart + ext_5prime
    return max(0, start0), end0, (1 if sstart <= send else -1)


def _tool_versions() -> dict:
    def _v(cmd: list[str], parse):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
            return parse((r.stdout or "") + (r.stderr or ""))
        except Exception:
            return "unknown"

    return {
        "minimap2":  _v(["minimap2", "--version"],  lambda s: s.strip().split()[0]),
        "samtools":  _v(["samtools", "--version"],  lambda s: s.split("\n")[0].split()[-1]),
        "tblastn":   _v(["tblastn",  "-version"],   lambda s: s.split("\n")[0].split()[-1]),
        "mafft":     _v(["mafft",    "--version"],  lambda s: s.strip().split("\n")[0]),
        "pysam":     pysam.__version__,
        "biopython": __import__("Bio").__version__,
    }


def _make_entry(
    gene: str,
    prob_type: str,
    action: str,
    versions: dict,
    cmds: list[str],
    stops_bef: int,
    stops_aft: int,
    evidence_extra: dict | None = None,
    candidate: dict | None = None,
) -> dict:
    return {
        "timestamp": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
        "gene": gene,
        "problem": prob_type,
        "evidence": {
            "type": "reads_tblastn_mafft_consensus",
            "stop_codons_before": stops_bef,
            "stop_codons_after":  stops_aft,
            **(evidence_extra or {}),
        },
        "candidate": candidate,
        "action": action,
        "tools": versions,
        "commands": cmds,
        "mitocurator_version": "0.1.0-dev",
    }


def _append_audit(audit_path: Path, entry: dict) -> None:
    with open(audit_path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
