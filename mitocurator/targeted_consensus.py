from __future__ import annotations
from pathlib import Path
from time import perf_counter
import statistics
import random
import gzip
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.Align import PairwiseAligner

from .utils import ensure_dir, safe_get, run_cmd, get_genetic_code
from .io import read_record
from .read_support import resolve_read_sets

TYPE_TO_PRESET = {
    "pacbio_hifi": ("map-hifi", "single_fastq"), "hifi": ("map-hifi", "single_fastq"),
    "pacbio_clr": ("map-pb", "single_fastq"), "clr": ("map-pb", "single_fastq"),
    "ont": ("map-ont", "single_fastq"), "nanopore": ("map-ont", "single_fastq"),
    "illumina_pe": ("sr", "paired_fastq"), "illumina": ("sr", "paired_fastq"), "pe": ("sr", "paired_fastq"),
    "illumina_se": ("sr", "single_fastq"), "se": ("sr", "single_fastq"),
}


def _read_fastq(path: Path):
    reads = []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        while True:
            h = fh.readline().rstrip()
            if not h:
                break
            s = fh.readline().rstrip(); fh.readline(); q = fh.readline().rstrip()
            reads.append((h, s, q))
    return reads


def _write_fastq(path: Path, reads):
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for h, s, q in reads:
            out.write(f"{h}\n{s}\n+\n{q}\n")


def _resolve_readset_type(config: dict, read_set: str) -> str:
    for rs in resolve_read_sets(config):
        if str(rs.get("name")) == read_set:
            rt = str(rs.get("type", "")).lower()
            if rt in TYPE_TO_PRESET:
                return rt
    for grp in (safe_get(config, ["reads", "long"], []) or []) + (safe_get(config, ["reads", "short"], []) or []):
        if isinstance(grp, dict) and str(grp.get("name")) == read_set:
            rt = str(grp.get("type", "")).lower()
            if rt in TYPE_TO_PRESET:
                return rt
    n = read_set.lower()
    if any(k in n for k in ["hifi", "pacbio"]): return "pacbio_hifi"
    if "clr" in n: return "pacbio_clr"
    if any(k in n for k in ["ont", "nanopore"]): return "ont"
    if "illumina" in n or "pe" in n: return "illumina_pe"
    return "illumina_se"


def _scan_orfs(seq: str, code: int, min_orf_nt: int, six_frames=True):
    cands = []
    strands = [("+", seq), ("-", str(Seq(seq).reverse_complement()))] if six_frames else [("+", seq)]
    stop_codons = {"TAA", "TAG", "AGA", "AGG"} if code == 5 else {"TAA", "TAG", "TGA"}
    oid = 1
    for strand, s in strands:
        for frame in (0, 1, 2):
            i = frame
            while i + 3 <= len(s):
                j = i
                while j + 3 <= len(s) and s[j:j+3] not in stop_codons:
                    j += 3
                nt = s[i:j]
                if len(nt) >= min_orf_nt:
                    aa = str(Seq(nt).translate(table=code, to_stop=False))
                    cands.append({"best_orf_id": f"orf{oid}", "orf_start": i+1, "orf_end": j, "orf_strand": strand, "orf_frame": frame, "orf_nt": nt,
                                  "orf_length_nt": len(nt), "orf_length_aa": len(aa), "internal_stop_count": aa[:-1].count("*"), "terminal_stop": "yes" if aa.endswith("*") else "no", "aa": aa})
                    oid += 1
                i = j + 3
    return cands


def _reference_protein(config: dict, gene: str, code: int):
    ref = safe_get(config, ["reference", "genbank"], None) or safe_get(config, ["input", "reference_gb"], None)
    if not ref or not Path(ref).exists():
        return None
    rec, _ = read_record(ref)
    for f in rec.features:
        if f.type != "CDS":
            continue
        g = "".join(f.qualifiers.get("gene", [""])).upper()
        if g == gene.upper():
            nt = str(f.extract(rec.seq)).upper()
            nt = nt[: (len(nt)//3)*3]
            return str(Seq(nt).translate(table=code, to_stop=False))
    return None


def _align_metrics(query_aa: str, ref_aa: str):
    if not query_aa or not ref_aa:
        return ".", ".", ".", "."
    al = PairwiseAligner(mode="global")
    aln = al.align(query_aa, ref_aa)[0]
    score = float(aln.score)
    q, r = aln[0], aln[1]
    matches = sum(1 for a, b in zip(q, r) if a == b and a != "-" and b != "-")
    aligned = sum(1 for a, b in zip(q, r) if a != "-" and b != "-")
    pid = round(100.0 * matches / max(aligned, 1), 2)
    covq = round(100.0 * aligned / max(len(query_aa), 1), 2)
    covr = round(100.0 * aligned / max(len(ref_aa), 1), 2)
    return pid, covq, covr, score


def _reference_protein(config: dict, gene: str, code: int):
    ref = safe_get(config, ["reference", "genbank"], None) or safe_get(config, ["input", "reference_gb"], None)
    if not ref or not Path(ref).exists():
        return None
    rec, _ = read_record(ref)
    for f in rec.features:
        if f.type != "CDS":
            continue
        g = "".join(f.qualifiers.get("gene", [""])).upper()
        if g == gene.upper():
            nt = str(f.extract(rec.seq)).upper()
            nt = nt[: (len(nt)//3)*3]
            return str(Seq(nt).translate(table=code, to_stop=False))
    return None


def _align_metrics(query_aa: str, ref_aa: str):
    if not query_aa or not ref_aa:
        return ".", ".", ".", "."
    al = PairwiseAligner(mode="global")
    aln = al.align(query_aa, ref_aa)[0]
    score = float(aln.score)
    q, r = aln[0], aln[1]
    matches = sum(1 for a, b in zip(q, r) if a == b and a != "-" and b != "-")
    aligned_q = sum(1 for a, b in zip(q, r) if a != "-" and b != "-")
    pid = round(100.0 * matches / max(aligned_q, 1), 2)
    covq = round(100.0 * aligned_q / max(len(query_aa), 1), 2)
    covr = round(100.0 * aligned_q / max(len(ref_aa), 1), 2)
    return pid, covq, covr, score


def run_targeted_consensus(config, root: Path, refinement_dir: Path, reconstruction_pools_dir: Path, outdir: Path):
    import pysam
    outdir = ensure_dir(outdir)
    cfa = ensure_dir(outdir / "consensus_fasta")
    tsv_in = reconstruction_pools_dir / "reconstruction_pools.tsv"
    tsv_out = outdir / "targeted_consensus.tsv"
    md_out = outdir / "targeted_consensus.md"
    cols = ["target_id","target_type","gene","read_set","pool_type","reference_used","consensus_fasta","consensus_length","n_bases","ambiguous_bases","ambiguous_fraction","mean_depth","min_depth","max_depth","best_orf_id","num_orfs_found","best_orf_selection_reason","orf_start","orf_end","orf_strand","orf_frame","orf_length_nt","orf_length_aa","internal_stop_count","terminal_stop","reference_gene","reference_aa_length","percent_identity","aligned_coverage_consensus","aligned_coverage_reference","alignment_score","recommendation","priority","reads_available","reads_used_for_consensus","downsampled","max_depth_per_position","elapsed_prepare_s","elapsed_map_s","elapsed_sort_s","elapsed_index_s","elapsed_pileup_s","elapsed_orf_s","elapsed_reference_compare_s","elapsed_write_s","elapsed_total_s","comment"]
    if not tsv_in.exists():
        pd.DataFrame(columns=cols).to_csv(tsv_out, sep="\t", index=False)
        md_out.write_text("# Targeted consensus\n\nNo reconstruction pools TSV found.\n", encoding="utf-8")
        return outdir

    df = pd.read_csv(tsv_in, sep="\t").fillna(".")
    pool_type = str(safe_get(config, ["targeted_consensus", "pool_type"], "combined"))
    df = df[df["pool_type"] == pool_type]
    tfilter = safe_get(config, ["targeted_consensus", "target_filter"], None)
    if tfilter:
        df = df[df["target_id"].astype(str).str.contains(str(tfilter), case=False) | df["gene"].astype(str).str.contains(str(tfilter), case=False)]
    rsfilter = safe_get(config, ["targeted_consensus", "read_set_filter"], None)
    if rsfilter:
        df = df[df["read_set"].astype(str).str.lower() == str(rsfilter).lower()]
    max_targets = safe_get(config, ["targeted_consensus", "max_targets"], None)
    if max_targets not in (None, "", "null"):
        keep = list(dict.fromkeys(df["target_id"].tolist()))[: int(max_targets)]
        df = df[df["target_id"].isin(keep)]

    min_bq = int(safe_get(config, ["targeted_consensus", "min_base_quality"], 20))
    min_mq = int(safe_get(config, ["targeted_consensus", "min_mapping_quality"], 20))
    min_depth = int(safe_get(config, ["targeted_consensus", "min_depth"], 5))
    maj = float(safe_get(config, ["targeted_consensus", "majority_threshold"], 0.7))
    minimap2_threads = int(safe_get(config, ["targeted_consensus", "minimap2_threads"], 4))
    samtools_threads = int(safe_get(config, ["targeted_consensus", "samtools_threads"], 4))
    max_reads = safe_get(config, ["targeted_consensus", "max_reads"], 10000)
    max_reads = None if max_reads in (None, "", "null") else int(max_reads)
    max_depth_pp = int(safe_get(config, ["targeted_consensus", "max_depth_per_position"], 1000))
    min_orf_nt = int(safe_get(config, ["targeted_consensus", "min_orf_nt"], 150))
    max_orf_aa_factor = float(safe_get(config, ["targeted_consensus", "max_orf_aa_factor"], 1.5))
    scan_six_frames = bool(safe_get(config, ["targeted_consensus", "scan_six_frames"], True))
    code = get_genetic_code(config, default=5)
    reuse = bool(safe_get(config, ["targeted_consensus", "reuse_existing_outputs"], True))
    rng = random.Random(int(safe_get(config, ["targeted_consensus", "random_seed"], 42)))

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
        t_total = perf_counter()
        tid, rs = str(r.target_id), str(r.read_set)
        if tid not in targets:
            continue
        rt = _resolve_readset_type(config, rs)
        preset, inferred_fmt = TYPE_TO_PRESET.get(rt, ("sr", "single_fastq"))
        task_fmt = str(r.output_format) if str(r.output_format) in {"single_fastq", "paired_fastq", "interleaved_fastq"} else inferred_fmt
        tdir = ensure_dir(outdir / tid / rs)
        adir = ensure_dir(tdir / "alignment")
        dsdir = ensure_dir(tdir / "reads_downsampled")
        local_ref = tdir / "local_ref.fa"
        _, s0, e0 = targets[tid]
        local_seq = str(rec.seq[s0:e0]).upper()
        SeqIO.write([SeqIO.SeqRecord(Seq(local_seq), id=tid, description="")], str(local_ref), "fasta")
        out_bam = adir / "pool_to_local_ref.bam"; out_bai = Path(str(out_bam) + ".bai")
        cpath = cfa / f"{tid}.{rs}.{pool_type}.consensus.fasta"

        # prepare reads
        t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=prepare_reads")
        fq, fq1, fq2 = str(r.output_fastq), str(r.output_fastq_r1), str(r.output_fastq_r2)
        use_fq, use_fq1, use_fq2 = fq, fq1, fq2
        reads_available = reads_used = 0; downsampled = "no"
        if max_reads:
            if task_fmt == "paired_fastq" and fq1 != "." and fq2 != ".":
                r1, r2 = _read_fastq(Path(fq1)), _read_fastq(Path(fq2)); n = min(len(r1), len(r2)); reads_available = n
                if n > max_reads:
                    idx = sorted(rng.sample(range(n), max_reads)); downsampled = "yes"
                    r1s = [r1[i] for i in idx]; r2s = [r2[i] for i in idx]
                    use_fq1 = str(dsdir / f"{tid}.{rs}.R1.ds.fastq.gz"); use_fq2 = str(dsdir / f"{tid}.{rs}.R2.ds.fastq.gz")
                    _write_fastq(Path(use_fq1), r1s); _write_fastq(Path(use_fq2), r2s); reads_used = len(r1s)
                else: reads_used = n
            else:
                rr = _read_fastq(Path(fq)); reads_available = len(rr)
                if len(rr) > max_reads:
                    idx = sorted(rng.sample(range(len(rr)), max_reads)); downsampled = "yes"
                    rrs = [rr[i] for i in idx]; use_fq = str(dsdir / f"{tid}.{rs}.ds.fastq.gz"); _write_fastq(Path(use_fq), rrs); reads_used = len(rrs)
                else: reads_used = len(rr)
        el_prep = perf_counter() - t0

        # map/sort/index
        el_map = el_sort = el_index = 0.0
        if not (reuse and out_bam.exists() and out_bai.exists() and cpath.exists()):
            t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=map")
            if task_fmt == "paired_fastq" and use_fq1 != "." and use_fq2 != ".":
                cmd = f"minimap2 -t {minimap2_threads} -ax {preset} {local_ref} {use_fq1} {use_fq2}"
            else:
                cmd = f"minimap2 -t {minimap2_threads} -ax {preset} {local_ref} {use_fq}"
            sam = adir / "tmp.sam"; run_cmd(["bash", "-lc", f"{cmd} > {sam}"], check=False); el_map = perf_counter() - t0
            t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=sort")
            run_cmd(["samtools", "sort", "-@", str(samtools_threads), "-o", str(out_bam), str(sam)], check=False); el_sort = perf_counter() - t0
            t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=index")
            run_cmd(["samtools", "index", str(out_bam)], check=False); el_index = perf_counter() - t0
            if sam.exists(): sam.unlink()
        else:
            print(f"[targeted_consensus] Reusing existing consensus for target={tid} read_set={rs}")

        # pileup
        t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=pileup")
        aln = pysam.AlignmentFile(str(out_bam), "rb") if out_bam.exists() else None
        cons = []; depths = []
        if aln is not None:
            for p in range(len(local_seq)):
                countsA=countsC=countsG=countsT=0; depth=0
                for col in aln.pileup(tid, p, p+1, truncate=True, stepper="all", max_depth=max_depth_pp):
                    if col.reference_pos != p: continue
                    for pr in col.pileups:
                        if pr.is_del or pr.is_refskip: continue
                        a=pr.alignment
                        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < min_mq: continue
                        qpos=pr.query_position
                        if qpos is None or qpos >= len(a.query_sequence): continue
                        bq = a.query_qualities[qpos] if a.query_qualities is not None else 40
                        if bq < min_bq: continue
                        b=a.query_sequence[qpos].upper(); depth += 1
                        if b=="A": countsA += 1
                        elif b=="C": countsC += 1
                        elif b=="G": countsG += 1
                        elif b=="T": countsT += 1
                depths.append(depth)
                if depth < min_depth: cons.append("N"); continue
                cts = {"A":countsA,"C":countsC,"G":countsG,"T":countsT}
                mb, mc = max(cts.items(), key=lambda x: x[1])
                cons.append(mb if (mc/depth) >= maj else "N")
            aln.close()
        cseq = "".join(cons); el_pile = perf_counter() - t0
        print(f"[targeted_consensus] target={tid} read_set={rs} step=pileup elapsed_s={round(el_pile,2)}")

        # orf
        t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=orf")
        orfs = _scan_orfs(cseq, code, min_orf_nt=min_orf_nt, six_frames=scan_six_frames)
        ref_aa = _reference_protein(config, str(r.gene), code)
        ref_len = len(ref_aa) if ref_aa else None
        filtered = orfs
        reason = "longest_min_stops"
        if ref_len:
            max_allowed = int(ref_len * max_orf_aa_factor)
            filt2 = [o for o in orfs if o["orf_length_aa"] <= max_allowed]
            if filt2:
                filtered = filt2; reason = "length_filtered_vs_reference"
        best = None; best_align = (-1, -1, -1, -1)
        if ref_aa and filtered:
            for o in filtered:
                pid,cq,cr,sc = _align_metrics(o["aa"], ref_aa)
                score = float(sc if sc != "." else -1)
                if score > best_align[3]:
                    best_align = (pid,cq,cr,score); best=o
            reason = "best_similarity_to_reference"
        elif filtered:
            best = sorted(filtered, key=lambda x: (-x["orf_length_aa"], x["internal_stop_count"]))[0]
        el_orf = perf_counter() - t0
        print(f"[targeted_consensus] target={tid} read_set={rs} step=orf elapsed_s={round(el_orf,2)}")

        # reference compare timing
        t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=reference_compare")
        pid=covc=covr=score=".";
        if best and ref_aa:
            pid,covc,covr,score = _align_metrics(best["aa"], ref_aa)
        el_ref = perf_counter() - t0
        print(f"[targeted_consensus] target={tid} read_set={rs} step=reference_compare elapsed_s={round(el_ref,2)}")

        # recommendation
        recm, pri = "MANUAL_REVIEW", "MEDIUM"
        if not best:
            recm = "NO_RELIABLE_CONSENSUS"
        elif statistics.mean(depths) if depths else 0 < min_depth:
            recm = "LOW_DEPTH_CONSENSUS"
        elif (cseq.count("N")/len(cseq)) if cseq else 1 > 0.3:
            recm = "CONSENSUS_AMBIGUOUS"
        elif str(r.target_type) == "missing_gene_candidate":
            if pid != "." and float(pid) >= 40 and float(covr) >= 50 and best["internal_stop_count"] == 0:
                recm, pri = "GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "HIGH"
            elif pid != "." and float(pid) >= 30 and float(covr) >= 25 and best["internal_stop_count"] == 0:
                recm, pri = "PARTIAL_GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "HIGH"
            elif pid != "." and best["internal_stop_count"] > 0:
                recm = "CONSENSUS_HAS_ORF_WITH_STOPS"
        elif "correction" in str(r.target_type) and best and best["internal_stop_count"] == 0:
            recm, pri = "LOCAL_CORRECTION_SUPPORTED_BY_CONSENSUS", "HIGH"
        elif "problematic" in str(r.target_type) and best and best["internal_stop_count"] > 0:
            recm, pri = "STOP_CONFIRMED_BY_CONSENSUS", "HIGH"

        t0 = perf_counter()
        SeqIO.write([SeqIO.SeqRecord(Seq(cseq), id=f"{tid}.{rs}", description="")], str(cfa / f"{tid}.{rs}.{pool_type}.consensus.fasta"), "fasta")
        el_write = perf_counter() - t0

        b = best or {"best_orf_id": ".", "orf_start": ".", "orf_end": ".", "orf_strand": ".", "orf_frame": ".", "orf_length_nt": 0, "orf_length_aa": 0, "internal_stop_count": ".", "terminal_stop": "."}
        rows.append({"target_id": tid, "target_type": r.target_type, "gene": r.gene, "read_set": rs, "pool_type": pool_type,
                     "reference_used": str(local_ref), "consensus_fasta": str(cfa / f"{tid}.{rs}.{pool_type}.consensus.fasta"),
                     "consensus_length": len(cseq), "n_bases": len(cseq)-cseq.count("N"), "ambiguous_bases": cseq.count("N"),
                     "ambiguous_fraction": round((cseq.count("N")/len(cseq)) if cseq else 1,4), "mean_depth": round(statistics.mean(depths),2) if depths else 0,
                     "min_depth": min(depths) if depths else 0, "max_depth": max(depths) if depths else 0,
                     "best_orf_id": b["best_orf_id"], "num_orfs_found": len(orfs), "best_orf_selection_reason": reason,
                     "orf_start": b["orf_start"], "orf_end": b["orf_end"], "orf_strand": b["orf_strand"], "orf_frame": b["orf_frame"],
                     "orf_length_nt": b["orf_length_nt"], "orf_length_aa": b["orf_length_aa"], "internal_stop_count": b["internal_stop_count"],
                     "terminal_stop": b["terminal_stop"], "reference_gene": str(r.gene), "reference_aa_length": ref_len if ref_len else ".",
                     "percent_identity": pid, "aligned_coverage_consensus": covc, "aligned_coverage_reference": covr, "alignment_score": score,
                     "recommendation": recm, "priority": pri, "reads_available": reads_available, "reads_used_for_consensus": reads_used,
                     "downsampled": downsampled, "max_depth_per_position": max_depth_pp,
                     "elapsed_prepare_s": round(el_prep,3), "elapsed_map_s": round(el_map,3), "elapsed_sort_s": round(el_sort,3), "elapsed_index_s": round(el_index,3),
                     "elapsed_pileup_s": round(el_pile,3), "elapsed_orf_s": round(el_orf,3), "elapsed_reference_compare_s": round(el_ref,3),
                     "elapsed_write_s": round(el_write,3), "elapsed_total_s": round(perf_counter()-t_total,3),
                     "comment": f"preset={preset}; type={rt}"})
        print(f"[targeted_consensus] target={tid} read_set={rs} step=done elapsed_s={round(perf_counter()-t_total,2)}")

    out_df = pd.DataFrame(rows, columns=cols)
    out_df.to_csv(tsv_out, sep="\t", index=False)
    with open(md_out, "w", encoding="utf-8") as md:
        md.write("# Targeted consensus\n\n")
        md.write(f"- total consensos: {len(out_df)}\n\n")
        sections = [
            ("GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "Candidatos de genes ausentes suportados"),
            ("PARTIAL_GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "Candidatos parciais suportados"),
            ("CONSENSUS_HAS_ORF_WITH_STOPS", "Consensos com ORF semelhante e stops"),
            ("MANUAL_REVIEW", "Casos para revisão manual"),
        ]
        for key, title in sections:
            md.write(f"## {title}\n\n")
            sub = out_df[out_df["recommendation"] == key] if not out_df.empty else pd.DataFrame()
            if sub.empty:
                md.write("- Nenhum\n\n")
            else:
                for rr in sub.itertuples():
                    md.write(f"- {rr.target_id} [{rr.read_set}] pid={rr.percent_identity} cov_ref={rr.aligned_coverage_reference}\n")
                md.write("\n")
        md.write("> Nenhum GenBank foi alterado.\n")
    return outdir
