from __future__ import annotations
from collections import Counter
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from .utils import ensure_dir, safe_get, run_cmd


def _load_pysam():
    try:
        import pysam
        return pysam
    except ImportError as exc:
        raise RuntimeError(
            "The read_support step requires the Python package 'pysam'. "
            "Install it with: mamba install -c bioconda -c conda-forge pysam"
        ) from exc
from .io import read_record
from .utils import get_genetic_code


def _extract_ref_fasta(refined_gb: Path, out_fa: Path):
    rec, _ = read_record(refined_gb)
    SeqIO.write(rec, str(out_fa), "fasta")


def _build_hifi_bam(config: dict, refined_fa: Path, outdir: Path):
    reads = safe_get(config, ["read_support", "hifi_reads"], []) or []
    if not reads:
        return None
    bam = outdir / "hifi_to_refined.bam"
    threads = str(int(safe_get(config, ["read_support", "threads"], 8)))
    cmd = f"minimap2 -t {threads} -ax map-hifi {refined_fa} {' '.join(reads)} | samtools sort -o {bam}"
    run_cmd(["bash", "-lc", cmd], check=False)
    run_cmd(["samtools", "index", str(bam)], check=False)
    return bam if bam.exists() else None


def _aa_for_codon(codon: str, code: int):
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "X"
    return str(Seq(codon).translate(table=code, to_stop=False))


def run_read_support(config: dict, refined_gb: Path, refinement_dir: Path, outdir: Path):
    pysam = _load_pysam()
    outdir = ensure_dir(outdir)
    stop_tsv = refinement_dir / "problematic_cds_stop_context.tsv"
    support_tsv = outdir / "problematic_stop_read_support.tsv"
    var_tsv = outdir / "problematic_stop_variants.tsv"
    summary_md = outdir / "read_support_summary.md"

    sup_cols = ["gene","seqid","cds_start","cds_end","strand","stop_aa_position","genomic_codon_start","genomic_codon_end","reference_codon","translation_table","read_depth_min","reads_covering_full_codon","major_read_codon","major_read_codon_count","major_read_codon_frequency","major_read_amino_acid","would_remove_stop","alternative_codons","alternative_codon_counts","recommendation","comment"]
    var_cols = ["gene","seqid","genomic_position","reference_base","depth","A","C","G","T","N","major_base","major_base_frequency","minor_bases","comment"]

    if not stop_tsv.exists():
        pd.DataFrame(columns=sup_cols).to_csv(support_tsv, sep="\t", index=False)
        pd.DataFrame(columns=var_cols).to_csv(var_tsv, sep="\t", index=False)
        summary_md.write_text("# Read support summary\n\nNo problematic stop context file found.\n", encoding="utf-8")
        return outdir

    refined_fa = outdir / "refined.fa"
    _extract_ref_fasta(refined_gb, refined_fa)
    bam = _build_hifi_bam(config, refined_fa, outdir) if bool(safe_get(config, ["read_support", "use_hifi"], True)) else None
    code = get_genetic_code(config, default=5)
    flank = int(safe_get(config, ["read_support", "flank_bp"], 20))

    stops = pd.read_csv(stop_tsv, sep="\t")
    support_rows, var_rows = [], []

    aln = pysam.AlignmentFile(str(bam), "rb") if bam and bam.exists() else None
    for _, r in stops.iterrows():
        g1, g2 = int(r["genomic_codon_start"]), int(r["genomic_codon_end"])
        pos = [g1-1, g1, g2-1] if g2 - g1 == 2 else [g1-1, g1, g1+1]
        ref_codon = str(r.get("codon", "NNN"))
        codon_counter = Counter()
        depth_per_pos = []

        if aln is not None:
            for p in [g1, g1+1, g1+2]:
                col = aln.pileup(str(r["seqid"]), p-1, p, truncate=True, stepper="all")
                counts = Counter()
                depth = 0
                for c in col:
                    if c.reference_pos != p-1:
                        continue
                    for pr in c.pileups:
                        if pr.is_del or pr.is_refskip:
                            continue
                        b = pr.alignment.query_sequence[pr.query_position].upper()
                        if b not in "ACGT":
                            b = "N"
                        counts[b] += 1
                        depth += 1
                depth_per_pos.append(depth)
                major = counts.most_common(1)[0][0] if counts else "N"
                freq = (counts[major] / depth) if depth else 0
                minors = ",".join([f"{k}:{v}" for k, v in counts.items() if k != major])
                var_rows.append({"gene": r["gene"], "seqid": r["seqid"], "genomic_position": p, "reference_base": ref_codon[p-g1] if len(ref_codon)==3 else "N", "depth": depth, "A": counts.get("A",0), "C": counts.get("C",0), "G": counts.get("G",0), "T": counts.get("T",0), "N": counts.get("N",0), "major_base": major, "major_base_frequency": round(freq,3), "minor_bases": minors if minors else ".", "comment": "pileup base counts"})

            for read in aln.fetch(str(r["seqid"]), max(0, g1-1-flank), g2+flank):
                qpos = []
                for rp, qp in read.get_aligned_pairs(matches_only=False):
                    if rp in (g1-1, g1, g1+1) and qp is not None:
                        qpos.append((rp, qp))
                if len(qpos) == 3:
                    qpos = sorted(qpos)
                    cod = ''.join(read.query_sequence[qp].upper() if read.query_sequence[qp].upper() in 'ACGT' else 'N' for _, qp in qpos)
                    codon_counter[cod] += 1

        full_cov = sum(codon_counter.values())
        major_codon, major_count = (codon_counter.most_common(1)[0] if codon_counter else ("NNN", 0))
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

        support_rows.append({
            "gene": r["gene"], "seqid": r["seqid"], "cds_start": r["cds_start"], "cds_end": r["cds_end"], "strand": r["strand"], "stop_aa_position": r["stop_aa_position"], "genomic_codon_start": g1, "genomic_codon_end": g2, "reference_codon": ref_codon, "translation_table": code,
            "read_depth_min": min(depth_per_pos) if depth_per_pos else 0, "reads_covering_full_codon": full_cov,
            "major_read_codon": major_codon, "major_read_codon_count": major_count, "major_read_codon_frequency": round(major_freq,3), "major_read_amino_acid": major_aa,
            "would_remove_stop": str(bool(would_remove_stop)).lower(),
            "alternative_codons": ";".join([k for k, _ in codon_counter.items() if k != major_codon]) if codon_counter else ".",
            "alternative_codon_counts": ";".join([f"{k}:{v}" for k, v in codon_counter.items() if k != major_codon]) if codon_counter else ".",
            "recommendation": rec,
            "comment": "read-level stop codon support (diagnostic-only)",
        })

    pd.DataFrame(support_rows, columns=sup_cols).to_csv(support_tsv, sep="\t", index=False)
    pd.DataFrame(var_rows, columns=var_cols).to_csv(var_tsv, sep="\t", index=False)

    sup_df = pd.DataFrame(support_rows)
    total = len(sup_df)
    supported = int((sup_df["recommendation"] == "STOP_SUPPORTED_BY_READS").sum()) if total else 0
    rescue = int((sup_df["recommendation"] == "POSSIBLE_POLISHING_OR_ASSEMBLY_ERROR").sum()) if total else 0
    pri = sup_df.sort_values(["gene", "major_read_codon_frequency"], ascending=[True, False]) if total else pd.DataFrame()
    with open(summary_md, "w", encoding="utf-8") as md:
        md.write("# Read support summary\n\n")
        md.write(f"- Total stops evaluated: {total}\n")
        md.write(f"- Stops supported by reads: {supported}\n")
        md.write(f"- Stops with possible rescue codon: {rescue}\n\n")
        md.write("## Prioritized by gene\n\n")
        if pri.empty:
            md.write("No read support rows generated.\n")
        else:
            for _, rr in pri.iterrows():
                md.write(f"- {rr['gene']} pos {rr['stop_aa_position']}: {rr['recommendation']} (major={rr['major_read_codon']} freq={rr['major_read_codon_frequency']})\n")
        md.write("\n> No automatic annotation change was applied.\n")

    return outdir
