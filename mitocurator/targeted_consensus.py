from __future__ import annotations
from pathlib import Path
import statistics
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

from .utils import ensure_dir, safe_get, run_cmd, get_genetic_code
from .io import read_record


def _translate_metrics(nt: str, code: int):
    nt = nt[: (len(nt) // 3) * 3]
    aa = str(Seq(nt).translate(table=code, to_stop=False)) if nt else ""
    return aa, aa[:-1].count("*") if aa else 0, "yes" if aa.endswith("*") else "no"


def _best_orf(seq: str, code: int):
    best = ("", 0, 0, "+", 0, 0, "no")
    for strand, s in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        for frame in (0, 1, 2):
            nt = s[frame:]
            nt = nt[: (len(nt)//3)*3]
            aa, internal, terminal = _translate_metrics(nt, code)
            cand = (nt, 1+frame, frame+len(nt), strand, frame, len(aa), internal, terminal)
            if len(aa) > best[5] or (len(aa)==best[5] and internal < best[6]):
                best = cand
    return best


def run_targeted_consensus(config, root: Path, refinement_dir: Path, reconstruction_pools_dir: Path, outdir: Path):
    import pysam  # lazy import
    outdir = ensure_dir(outdir)
    cfa = ensure_dir(outdir / "consensus_fasta")
    tsv_in = reconstruction_pools_dir / "reconstruction_pools.tsv"
    tsv_out = outdir / "targeted_consensus.tsv"
    md_out = outdir / "targeted_consensus.md"
    cols = ["target_id","target_type","gene","read_set","pool_type","reference_used","consensus_fasta","consensus_length","n_bases","ambiguous_bases","ambiguous_fraction","mean_depth","min_depth","max_depth","orf_start","orf_end","orf_strand","orf_frame","orf_length_nt","orf_length_aa","internal_stop_count","terminal_stop","reference_gene","reference_aa_length","percent_identity","aligned_coverage_consensus","aligned_coverage_reference","recommendation","priority","comment"]
    if not tsv_in.exists():
        pd.DataFrame(columns=cols).to_csv(tsv_out, sep="\t", index=False)
        md_out.write_text("# Targeted consensus\n\nNo reconstruction pools TSV found.\n", encoding="utf-8")
        return outdir

    df = pd.read_csv(tsv_in, sep="\t").fillna(".")
    pool_type = str(safe_get(config, ["targeted_consensus", "pool_type"], "combined"))
    df = df[df["pool_type"] == pool_type]
    min_bq = int(safe_get(config, ["targeted_consensus", "min_base_quality"], 20))
    min_mq = int(safe_get(config, ["targeted_consensus", "min_mapping_quality"], 20))
    min_depth = int(safe_get(config, ["targeted_consensus", "min_depth"], 5))
    maj = float(safe_get(config, ["targeted_consensus", "majority_threshold"], 0.7))
    code = get_genetic_code(config, default=5)
    rec, _ = read_record(root / "05_refinement" / "refined.gb")
    targets = {}
    bed = root / "08_targeted_extraction" / "targets.bed"
    if bed.exists():
        with open(bed, "r", encoding="utf-8") as f:
            next(f, None)
            for ln in f:
                seqid, s0, e0, tid, *_ = ln.rstrip("\n").split("\t")
                targets[tid] = (seqid, int(s0), int(e0))

    rows = []
    for r in df.itertuples():
        tid = str(r.target_id); rs = str(r.read_set)
        tdir = ensure_dir(outdir / tid / rs)
        adir = ensure_dir(tdir / "alignment")
        if tid not in targets:
            continue
        seqid, s0, e0 = targets[tid]
        local_ref = tdir / "local_ref.fa"
        local_seq = str(rec.seq[s0:e0]).upper()
        SeqIO.write([SeqIO.SeqRecord(Seq(local_seq), id=tid, description="")], str(local_ref), "fasta")

        out_bam = adir / "pool_to_local_ref.bam"
        fq = str(r.output_fastq)
        fq1, fq2 = str(r.output_fastq_r1), str(r.output_fastq_r2)
        if str(r.output_format) == "paired_fastq" and fq1 != "." and fq2 != ".":
            cmd = f"minimap2 -ax sr {local_ref} {fq1} {fq2} | samtools sort -o {out_bam}"
        else:
            cmd = f"minimap2 -ax sr {local_ref} {fq} | samtools sort -o {out_bam}"
        run_cmd(["bash", "-lc", cmd], check=False)
        run_cmd(["samtools", "index", str(out_bam)], check=False)

        aln = pysam.AlignmentFile(str(out_bam), "rb") if out_bam.exists() else None
        consensus = []
        depths = []
        if aln is not None:
            for p in range(len(local_seq)):
                counts = {"A":0,"C":0,"G":0,"T":0}
                depth = 0
                for col in aln.pileup(tid, p, p+1, truncate=True, stepper="all"):
                    if col.reference_pos != p:
                        continue
                    for pr in col.pileups:
                        if pr.is_del or pr.is_refskip:
                            continue
                        a = pr.alignment
                        if a.mapping_quality < min_mq:
                            continue
                        qpos = pr.query_position
                        if qpos is None or qpos >= len(a.query_sequence):
                            continue
                        bq = a.query_qualities[qpos] if a.query_qualities is not None else 40
                        if bq < min_bq:
                            continue
                        b = a.query_sequence[qpos].upper()
                        if b in counts:
                            counts[b] += 1
                            depth += 1
                depths.append(depth)
                if depth < min_depth:
                    consensus.append("N"); continue
                mb, mc = max(counts.items(), key=lambda x: x[1])
                consensus.append(mb if (mc/depth) >= maj else "N")
            aln.close()
        cseq = "".join(consensus) if consensus else ""
        cpath = cfa / f"{tid}.{rs}.{pool_type}.consensus.fasta"
        SeqIO.write([SeqIO.SeqRecord(Seq(cseq), id=f"{tid}.{rs}", description="")], str(cpath), "fasta")

        nt, orf_s, orf_e, orf_strand, orf_frame, orf_aa_len, istops, terminal = _best_orf(cseq, code)
        ncount = cseq.count("N")
        amb_frac = (ncount / len(cseq)) if cseq else 1.0
        dmean = statistics.mean(depths) if depths else 0
        dmin = min(depths) if depths else 0
        dmax = max(depths) if depths else 0

        pid, covc, covr, ref_gene, ref_aa_len = ".", ".", ".", str(r.gene), "."
        if bool(safe_get(config, ["targeted_consensus", "compare_to_reference"], True)) and orf_aa_len > 0:
            # lightweight self-comparison placeholder using problematic_cds_reference_check when available
            ref_tsv = refinement_dir / "problematic_cds_reference_check.tsv"
            if ref_tsv.exists():
                rf = pd.read_csv(ref_tsv, sep="\t")
                hit = rf[rf["gene"].astype(str) == str(r.gene)]
                if not hit.empty:
                    ref_aa_len = int(hit.iloc[0].get("reference_aa_length", 0) or 0)
                    if ref_aa_len > 0:
                        pid = float(hit.iloc[0].get("percent_identity", 0) or 0)
                        covr = float(hit.iloc[0].get("aligned_coverage_reference", 0) or 0)
                        covc = min(100.0, round((orf_aa_len / ref_aa_len) * 100.0, 2))

        if dmean < min_depth:
            recm, pri = "LOW_DEPTH_CONSENSUS", "MEDIUM"
        elif amb_frac > 0.3:
            recm, pri = "CONSENSUS_AMBIGUOUS", "MEDIUM"
        elif str(r.target_type) == "missing_gene_candidate" and istops == 0 and orf_aa_len > 80:
            recm, pri = "GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "HIGH"
        elif "correction" in str(r.target_type) and istops == 0:
            recm, pri = "LOCAL_CORRECTION_SUPPORTED_BY_CONSENSUS", "HIGH"
        elif "problematic" in str(r.target_type) and istops > 0 and dmean >= min_depth:
            recm, pri = "STOP_CONFIRMED_BY_CONSENSUS", "HIGH"
        elif len(cseq) == 0:
            recm, pri = "NO_RELIABLE_CONSENSUS", "MEDIUM"
        else:
            recm, pri = "MANUAL_REVIEW", "MEDIUM"

        rows.append({"target_id": tid, "target_type": r.target_type, "gene": r.gene, "read_set": rs, "pool_type": pool_type, "reference_used": str(local_ref), "consensus_fasta": str(cpath), "consensus_length": len(cseq), "n_bases": len(cseq)-ncount, "ambiguous_bases": ncount, "ambiguous_fraction": round(amb_frac,4), "mean_depth": round(dmean,2), "min_depth": dmin, "max_depth": dmax, "orf_start": orf_s, "orf_end": orf_e, "orf_strand": orf_strand, "orf_frame": orf_frame, "orf_length_nt": len(nt), "orf_length_aa": orf_aa_len, "internal_stop_count": istops, "terminal_stop": terminal, "reference_gene": ref_gene, "reference_aa_length": ref_aa_len, "percent_identity": pid, "aligned_coverage_consensus": covc, "aligned_coverage_reference": covr, "recommendation": recm, "priority": pri, "comment": "diagnostic-only consensus"})

    out_df = pd.DataFrame(rows, columns=cols)
    out_df.to_csv(tsv_out, sep="\t", index=False)
    with open(md_out, "w", encoding="utf-8") as md:
        md.write("# Targeted consensus\n\n")
        md.write(f"- total consensos: {len(out_df)}\n\n")
        for k in ["GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "LOCAL_CORRECTION_SUPPORTED_BY_CONSENSUS", "STOP_CONFIRMED_BY_CONSENSUS", "CONSENSUS_AMBIGUOUS", "LOW_DEPTH_CONSENSUS", "NO_RELIABLE_CONSENSUS"]:
            md.write(f"- {k}: {int((out_df['recommendation']==k).sum()) if not out_df.empty else 0}\n")
        md.write("\n> Nenhum GenBank foi alterado.\n")
    return outdir
