from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .utils import ensure_dir, safe_get, run_cmd, get_genetic_code
from .io import read_record

DEFAULT_PRESET = {
    "pacbio_hifi": "map-hifi", "hifi": "map-hifi",
    "pacbio_clr": "map-pb", "clr": "map-pb",
    "ont": "map-ont", "nanopore": "map-ont",
    "illumina_pe": "sr", "illumina_se": "sr",
}


def resolve_read_sets(config):
    rs = safe_get(config, ["read_support", "read_sets"], None)
    if rs:
        out = []
        for i, r in enumerate(rs):
            name = r.get("name") or f"readset{i+1}"
            rtype = (r.get("type") or "pacbio_hifi").lower()
            mapper = r.get("mapper", "minimap2")
            preset = r.get("preset") or DEFAULT_PRESET.get(rtype, "map-hifi")
            reads = r.get("reads")
            if not reads:
                continue
            out.append({"name": name, "type": rtype, "mapper": mapper, "preset": preset, "reads": reads})
        return out

    # backward compatibility
    out = []
    if bool(safe_get(config, ["read_support", "use_hifi"], False)):
        hifi = safe_get(config, ["read_support", "hifi_reads"], []) or []
        if hifi:
            out.append({"name": "hifi", "type": "pacbio_hifi", "mapper": "minimap2", "preset": "map-hifi", "reads": hifi})
    if bool(safe_get(config, ["read_support", "use_illumina"], False)):
        il = safe_get(config, ["read_support", "illumina_reads"], {}) or {}
        if il.get("r1") and il.get("r2"):
            out.append({"name": "illumina", "type": "illumina_pe", "mapper": "minimap2", "preset": "sr", "reads": {"r1": il["r1"], "r2": il["r2"]}})
    return out


def _extract_ref_fasta(refined_gb: Path, out_fa: Path):
    rec, _ = read_record(refined_gb)
    SeqIO.write(rec, str(out_fa), "fasta")


def _build_bam_for_read_set(read_set: dict, refined_fa: Path, outdir: Path, threads: int):
    name, preset = read_set["name"], read_set["preset"]
    bam = outdir / f"{name}_to_refined.bam"
    reads = read_set["reads"]
    if isinstance(reads, dict):
        read_args = f"{reads.get('r1','')} {reads.get('r2','')}".strip()
    else:
        read_args = " ".join(reads)
    cmd = f"minimap2 -t {threads} -ax {preset} {refined_fa} {read_args} | samtools sort -o {bam}"
    run_cmd(["bash", "-lc", cmd], check=False)
    run_cmd(["samtools", "index", str(bam)], check=False)
    return bam if bam.exists() else None


def _aa_for_codon(codon: str, code: int):
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "X"
    return str(Seq(codon).translate(table=code, to_stop=False))


def run_read_support(config: dict, refined_gb: Path, refinement_dir: Path, outdir: Path):
    import pysam  # lazy import

    outdir = ensure_dir(outdir)
    stop_tsv = refinement_dir / "problematic_cds_stop_context.tsv"
    support_tsv = outdir / "problematic_stop_read_support.tsv"
    var_tsv = outdir / "problematic_stop_variants.tsv"
    summary_md = outdir / "read_support_summary.md"

    sup_cols = ["read_set","gene","seqid","cds_start","cds_end","strand","stop_aa_position","genomic_codon_start","genomic_codon_end","reference_codon","translation_table","read_depth_min","reads_covering_full_codon","major_read_codon","major_read_codon_count","major_read_codon_frequency","major_read_amino_acid","would_remove_stop","alternative_codons","alternative_codon_counts","recommendation","comment"]
    var_cols = ["read_set","gene","seqid","genomic_position","reference_base","depth","A","C","G","T","N","major_base","major_base_frequency","minor_bases","comment"]

    if not stop_tsv.exists():
        pd.DataFrame(columns=sup_cols).to_csv(support_tsv, sep="\t", index=False)
        pd.DataFrame(columns=var_cols).to_csv(var_tsv, sep="\t", index=False)
        summary_md.write_text("# Read support summary\n\nNo problematic stop context file found.\n", encoding="utf-8")
        return outdir

    read_sets = resolve_read_sets(config)
    refined_fa = outdir / "refined.fa"
    _extract_ref_fasta(refined_gb, refined_fa)
    threads = int(safe_get(config, ["read_support", "threads"], 8))
    bams = {rs["name"]: _build_bam_for_read_set(rs, refined_fa, outdir, threads) for rs in read_sets}

    code = get_genetic_code(config, default=5)
    flank = int(safe_get(config, ["read_support", "flank_bp"], 20))
    stops = pd.read_csv(stop_tsv, sep="\t")
    support_rows, var_rows = [], []

    for rs in read_sets:
        bam = bams.get(rs["name"])
        aln = pysam.AlignmentFile(str(bam), "rb") if bam and bam.exists() else None
        for _, r in stops.iterrows():
            g1, g2 = int(r["genomic_codon_start"]), int(r["genomic_codon_end"])
            ref_codon = str(r.get("codon", "NNN"))
            codon_counter, depth_per_pos = Counter(), []

            if aln is not None:
                for p in [g1, g1 + 1, g1 + 2]:
                    counts, depth = Counter(), 0
                    for c in aln.pileup(str(r["seqid"]), p - 1, p, truncate=True, stepper="all"):
                        if c.reference_pos != p - 1:
                            continue
                        for pr in c.pileups:
                            if pr.is_del or pr.is_refskip:
                                continue
                            b = pr.alignment.query_sequence[pr.query_position].upper()
                            b = b if b in "ACGT" else "N"
                            counts[b] += 1
                            depth += 1
                    depth_per_pos.append(depth)
                    major = counts.most_common(1)[0][0] if counts else "N"
                    freq = (counts[major] / depth) if depth else 0
                    minors = ",".join([f"{k}:{v}" for k, v in counts.items() if k != major])
                    var_rows.append({"read_set": rs["name"], "gene": r["gene"], "seqid": r["seqid"], "genomic_position": p, "reference_base": ref_codon[p-g1] if len(ref_codon)==3 else "N", "depth": depth, "A": counts.get("A",0), "C": counts.get("C",0), "G": counts.get("G",0), "T": counts.get("T",0), "N": counts.get("N",0), "major_base": major, "major_base_frequency": round(freq,3), "minor_bases": minors if minors else ".", "comment": "pileup base counts"})

                for read in aln.fetch(str(r["seqid"]), max(0, g1 - 1 - flank), g2 + flank):
                    qpos = [(rp, qp) for rp, qp in read.get_aligned_pairs(matches_only=False) if rp in (g1-1, g1, g1+1) and qp is not None]
                    if len(qpos) == 3:
                        qpos = sorted(qpos)
                        cod = ''.join(read.query_sequence[qp].upper() if read.query_sequence[qp].upper() in 'ACGT' else 'N' for _, qp in qpos)
                        codon_counter[cod] += 1

            full_cov = sum(codon_counter.values())
            major_codon, major_count = codon_counter.most_common(1)[0] if codon_counter else ("NNN", 0)
            major_freq = (major_count / full_cov) if full_cov else 0.0
            major_aa = _aa_for_codon(major_codon, code)
            would_remove_stop = (major_aa != "*") and major_codon != "NNN"

            if full_cov < 5:
                rec = "LOW_COVERAGE_REVIEW"
            elif major_codon == ref_codon and major_freq >= 0.8:
                rec = "STOP_SUPPORTED_BY_READS"
            elif major_aa != "*" and major_freq >= 0.7:
                rec = "POSSIBLE_POLISHING_OR_ASSEMBLY_ERROR"
            elif major_freq < 0.7:
                rec = "AMBIGUOUS_READ_SUPPORT"
            else:
                rec = "MANUAL_REVIEW"

            support_rows.append({"read_set": rs["name"], "gene": r["gene"], "seqid": r["seqid"], "cds_start": r["cds_start"], "cds_end": r["cds_end"], "strand": r["strand"], "stop_aa_position": r["stop_aa_position"], "genomic_codon_start": g1, "genomic_codon_end": g2, "reference_codon": ref_codon, "translation_table": code, "read_depth_min": min(depth_per_pos) if depth_per_pos else 0, "reads_covering_full_codon": full_cov, "major_read_codon": major_codon, "major_read_codon_count": major_count, "major_read_codon_frequency": round(major_freq,3), "major_read_amino_acid": major_aa, "would_remove_stop": str(bool(would_remove_stop)).lower(), "alternative_codons": ";".join([k for k, _ in codon_counter.items() if k != major_codon]) if codon_counter else ".", "alternative_codon_counts": ";".join([f"{k}:{v}" for k, v in codon_counter.items() if k != major_codon]) if codon_counter else ".", "recommendation": rec, "comment": "read-level stop codon support (diagnostic-only)"})

    pd.DataFrame(support_rows, columns=sup_cols).to_csv(support_tsv, sep="\t", index=False)
    pd.DataFrame(var_rows, columns=var_cols).to_csv(var_tsv, sep="\t", index=False)

    sup_df = pd.DataFrame(support_rows)
    with open(summary_md, "w", encoding="utf-8") as md:
        md.write("# Read support summary\n\n")
        if sup_df.empty:
            md.write("No read support rows generated.\n")
        else:
            for rs, g in sup_df.groupby("read_set"):
                md.write(f"## {rs}\n\n")
                md.write(f"- total de stops avaliados: {len(g)}\n")
                for k in ["STOP_SUPPORTED_BY_READS", "POSSIBLE_POLISHING_OR_ASSEMBLY_ERROR", "LOW_COVERAGE_REVIEW", "AMBIGUOUS_READ_SUPPORT"]:
                    md.write(f"- {k}: {int((g['recommendation']==k).sum())}\n")
                md.write("\n")
        md.write("> No automatic annotation change was applied.\n")
    return outdir
