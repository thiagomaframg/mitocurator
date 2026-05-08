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


def _parse_mpileup_bases(ref_base: str, bases: str):
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    i = 0
    ref_base = ref_base.upper()
    while i < len(bases):
        c = bases[i]
        if c == "^":
            i += 2; continue
        if c == "$":
            i += 1; continue
        if c in "+-":
            i += 1
            n = ""
            while i < len(bases) and bases[i].isdigit():
                n += bases[i]; i += 1
            if n: i += int(n)
            continue
        b = ref_base if c in "., " else c.upper()
        if b in counts: counts[b] += 1
        i += 1
    return counts


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
    cols = ["target_id","target_type","gene","read_set","pool_type","consensus_backend","reference_used","consensus_fasta","consensus_length","n_bases","ambiguous_bases","ambiguous_fraction","covered_bases","covered_fraction","low_depth_bases","low_depth_fraction","mean_depth","min_depth","max_depth","best_orf_id","num_orfs_found","best_orf_selection_reason","orf_start","orf_end","orf_strand","orf_frame","orf_length_nt","orf_length_aa","internal_stop_count","terminal_stop","reference_gene","reference_aa_length","percent_identity","aligned_coverage_consensus","aligned_coverage_reference","alignment_score","recommendation","priority","reads_available","reads_used_for_consensus","downsampled","max_depth_per_position","elapsed_prepare_s","elapsed_map_s","elapsed_sort_s","elapsed_index_s","elapsed_filter_s","elapsed_pileup_s","elapsed_orf_s","elapsed_reference_compare_s","elapsed_write_s","elapsed_total_s","comment"]
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
    consensus_backend = str(safe_get(config, ["targeted_consensus", "consensus_backend"], "count_coverage")).lower()
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
            run_cmd(["samtools", "index", "-@", str(samtools_threads), str(out_bam)], check=False); el_index = perf_counter() - t0
            if sam.exists(): sam.unlink()
        else:
            print(f"[targeted_consensus] Reusing existing consensus for target={tid} read_set={rs}")

        # filter
        t0 = perf_counter()
        filtered_bam = adir / "pool_to_local_ref.filtered.bam"
        run_cmd(["samtools", "view", "-@", str(samtools_threads), "-b", "-q", str(min_mq), "-F", "2308", str(out_bam), "-o", str(filtered_bam)], check=False)
        run_cmd(["samtools", "index", "-@", str(samtools_threads), str(filtered_bam)], check=False)
        el_filter = perf_counter() - t0

        # pileup
        t0 = perf_counter(); print(f"[targeted_consensus] target={tid} read_set={rs} step=pileup")
        aln = pysam.AlignmentFile(str(filtered_bam), "rb") if filtered_bam.exists() else None
        cons = []; depths = []
        if consensus_backend == "count_coverage" and aln is not None:
            A, C, G, T = aln.count_coverage(tid, start=0, stop=len(local_seq), quality_threshold=min_bq)
            for i in range(len(local_seq)):
                a, c, g, t = int(A[i]), int(C[i]), int(G[i]), int(T[i])
                depth = a + c + g + t
                depths.append(depth)
                if depth < min_depth:
                    cons.append("N"); continue
                cts = {"A": a, "C": c, "G": g, "T": t}
                mb, mc = max(cts.items(), key=lambda x: x[1])
                cons.append(mb if (mc / max(depth, 1)) >= maj else "N")
            aln.close()
        elif consensus_backend == "samtools" and filtered_bam.exists():
            mp = adir / "pileup.txt"
            print("[targeted_consensus] samtools mpileup does not support threads in this environment")
            run_cmd(["bash", "-lc", f"samtools mpileup -aa -Q {min_bq} -q {min_mq} -d {max_depth_pp} -f {local_ref} {filtered_bam} > {mp}"], check=False)
            pos = {}
            if mp.exists():
                with open(mp, "r", encoding="utf-8") as fh:
                    for ln in fh:
                        c = ln.rstrip("\n").split("\t")
                        if len(c) < 5:
                            continue
                        p = int(c[1]) - 1
                        depth = int(c[3]) if c[3].isdigit() else 0
                        counts = _parse_mpileup_bases(c[2], c[4])
                        pos[p] = (depth, counts)
            for p, rb in enumerate(local_seq):
                depth, cts = pos.get(p, (0, {"A": 0, "C": 0, "G": 0, "T": 0}))
                depths.append(depth)
                if depth < min_depth:
                    cons.append("N"); continue
                mb, mc = max(cts.items(), key=lambda x: x[1])
                cons.append(mb if (mc/max(depth,1)) >= maj else "N")
        elif aln is not None:
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
        mean_depth_obs = statistics.mean(depths) if depths else 0
        min_depth_obs = min(depths) if depths else 0
        n_bases = len(cseq) - cseq.count("N")
        amb_frac = (cseq.count("N")/len(cseq)) if cseq else 1
        covered_bases = sum(1 for d in depths if d >= 1)
        covered_fraction = (covered_bases / len(cseq)) if cseq else 0
        low_depth_bases = sum(1 for d in depths if d < min_depth)
        low_depth_fraction = (low_depth_bases / len(cseq)) if cseq else 1

        if not best:
            recm = "NO_RELIABLE_CONSENSUS"
        elif (mean_depth_obs < min_depth) or (n_bases == 0) or (covered_fraction < 0.80):
            recm = "LOW_DEPTH_CONSENSUS"
        elif amb_frac > 0.20:
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
        comment_extra = "max_depth_per_position not applied by count_coverage backend" if consensus_backend == "count_coverage" else ""
        note_zero = "; some positions have zero coverage" if (min_depth_obs == 0 and covered_fraction >= 0.80 and amb_frac <= 0.20) else ""
        rows.append({"target_id": tid, "target_type": r.target_type, "gene": r.gene, "read_set": rs, "pool_type": pool_type, "consensus_backend": consensus_backend,
                     "reference_used": str(local_ref), "consensus_fasta": str(cfa / f"{tid}.{rs}.{pool_type}.consensus.fasta"),
                     "consensus_length": len(cseq), "n_bases": n_bases, "ambiguous_bases": cseq.count("N"),
                     "ambiguous_fraction": round(amb_frac,4), "covered_bases": covered_bases, "covered_fraction": round(covered_fraction,4), "low_depth_bases": low_depth_bases, "low_depth_fraction": round(low_depth_fraction,4), "mean_depth": round(mean_depth_obs,2) if depths else 0,
                     "min_depth": min_depth_obs, "max_depth": max(depths) if depths else 0,
                     "best_orf_id": b["best_orf_id"], "num_orfs_found": len(orfs), "best_orf_selection_reason": reason,
                     "orf_start": b["orf_start"], "orf_end": b["orf_end"], "orf_strand": b["orf_strand"], "orf_frame": b["orf_frame"],
                     "orf_length_nt": b["orf_length_nt"], "orf_length_aa": b["orf_length_aa"], "internal_stop_count": b["internal_stop_count"],
                     "terminal_stop": b["terminal_stop"], "reference_gene": str(r.gene), "reference_aa_length": ref_len if ref_len else ".",
                     "percent_identity": pid, "aligned_coverage_consensus": covc, "aligned_coverage_reference": covr, "alignment_score": score,
                     "recommendation": recm, "priority": pri, "reads_available": reads_available, "reads_used_for_consensus": reads_used,
                     "downsampled": downsampled, "max_depth_per_position": max_depth_pp,
                     "elapsed_prepare_s": round(el_prep,3), "elapsed_map_s": round(el_map,3), "elapsed_sort_s": round(el_sort,3), "elapsed_index_s": round(el_index,3), "elapsed_filter_s": round(el_filter,3),
                     "elapsed_pileup_s": round(el_pile,3), "elapsed_orf_s": round(el_orf,3), "elapsed_reference_compare_s": round(el_ref,3),
                     "elapsed_write_s": round(el_write,3), "elapsed_total_s": round(perf_counter()-t_total,3),
                     "comment": f"preset={preset}; type={rt}" + (f"; {comment_extra}" if comment_extra else "") + note_zero})
        print(f"[targeted_consensus] target={tid} read_set={rs} step=done elapsed_s={round(perf_counter()-t_total,2)}")

    out_df = pd.DataFrame(rows, columns=cols)
    out_df.to_csv(tsv_out, sep="\t", index=False)
    # ranking for missing gene candidates
    rank_tsv = outdir / "best_missing_gene_candidates.tsv"
    rank_md = outdir / "best_missing_gene_candidates.md"
    sub = out_df[out_df["target_type"] == "missing_gene_candidate"].copy() if not out_df.empty else pd.DataFrame()
    if sub.empty:
        pd.DataFrame(columns=["gene","read_set","rank","target_id","recommendation","rank_score","selection_reason","orf_length_aa","reference_aa_length","internal_stop_count","percent_identity","aligned_coverage_reference","aligned_coverage_consensus","ambiguous_fraction","covered_fraction","mean_depth","consensus_fasta"]).to_csv(rank_tsv, sep="\t", index=False)
        rank_md.write_text("# Best missing-gene candidates\n\nNo missing-gene candidates.\n", encoding="utf-8")
    else:
        pri = {"GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS":7,"PARTIAL_GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS":6,"CONSENSUS_HAS_ORF_WITH_STOPS":5,"MANUAL_REVIEW":4,"CONSENSUS_AMBIGUOUS":3,"LOW_DEPTH_CONSENSUS":2,"NO_RELIABLE_CONSENSUS":1}
        for c in ["internal_stop_count","aligned_coverage_reference","percent_identity","ambiguous_fraction","covered_fraction","mean_depth","orf_length_aa","reference_aa_length"]:
            sub[c] = pd.to_numeric(sub[c], errors="coerce")
        sub["recommendation_priority"] = sub["recommendation"].map(pri).fillna(0)
        sub["len_delta"] = (sub["orf_length_aa"] - sub["reference_aa_length"]).abs()
        sub["rank_score"] = sub["recommendation_priority"]*10000 - sub["internal_stop_count"].fillna(999)*200 + sub["aligned_coverage_reference"].fillna(0)*10 + sub["percent_identity"].fillna(0)*5 - sub["len_delta"].fillna(999) - sub["ambiguous_fraction"].fillna(1)*100 + sub["covered_fraction"].fillna(0)*50 + sub["mean_depth"].fillna(0)*0.01
        out_rows = []
        for (g, rs), grp in sub.groupby(["gene","read_set"]):
            grp = grp.sort_values(["rank_score"], ascending=False).reset_index(drop=True)
            for i, rr in grp.iterrows():
                out_rows.append({"gene": g, "read_set": rs, "rank": i+1, "target_id": rr["target_id"], "recommendation": rr["recommendation"], "rank_score": round(rr["rank_score"],3), "selection_reason": "priority+stops+identity+coverage+length+ambiguity", "orf_length_aa": rr["orf_length_aa"], "reference_aa_length": rr["reference_aa_length"], "internal_stop_count": rr["internal_stop_count"], "percent_identity": rr["percent_identity"], "aligned_coverage_reference": rr["aligned_coverage_reference"], "aligned_coverage_consensus": rr["aligned_coverage_consensus"], "ambiguous_fraction": rr["ambiguous_fraction"], "covered_fraction": rr["covered_fraction"], "mean_depth": rr["mean_depth"], "consensus_fasta": rr["consensus_fasta"]})
        rank_df = pd.DataFrame(out_rows)
        rank_df.to_csv(rank_tsv, sep="\t", index=False)
        with open(rank_md, "w", encoding="utf-8") as md:
            md.write("# Best missing-gene candidates\n\n")
            for (g, rs), grp in rank_df.groupby(["gene","read_set"]):
                md.write(f"## {g} / {rs}\n\n")
                for rr in grp.itertuples():
                    md.write(f"- rank {rr.rank}: {rr.target_id} | {rr.recommendation} | pid={rr.percent_identity} cov_ref={rr.aligned_coverage_reference}\n")
                md.write("\n")
            md.write("> Nenhum GenBank foi alterado.\n")
    with open(md_out, "w", encoding="utf-8") as md:
        md.write("# Targeted consensus\n\n")
        md.write(f"- total consensos: {len(out_df)}\n\n")
        sections = [
            ("GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "Missing-gene candidates supported by consensus"),
            ("PARTIAL_GENE_CANDIDATE_SUPPORTED_BY_CONSENSUS", "Partial gene candidates supported by consensus"),
            ("CONSENSUS_HAS_ORF_WITH_STOPS", "Consensus with ORF similar to reference but with stops"),
            ("MANUAL_REVIEW", "Manual review cases"),
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
        md.write("## Best missing-gene candidates\n\n")
        md.write(f"- TSV: {rank_tsv}\n")
        md.write(f"- Markdown: {rank_md}\n\n")
        cross_tsv, cross_md = write_cross_readset_missing_gene_candidates(out_df, outdir, rank_tsv)
        md.write("## Cross-readset missing-gene ranking\n\n")
        md.write(f"- TSV: {cross_tsv}\n")
        md.write(f"- Markdown: {cross_md}\n\n")
        md.write("> Nenhum GenBank foi alterado.\n")
    return outdir


def write_cross_readset_missing_gene_candidates(targeted_df: pd.DataFrame, outdir: Path, best_rank_tsv: Path):
    cross_tsv = outdir / "cross_readset_missing_gene_candidates.tsv"
    cross_md = outdir / "cross_readset_missing_gene_candidates.md"
    cols = ["gene","target_id","combined_rank","combined_recommendation","read_sets_supporting","n_read_sets","best_rank","mean_rank","mean_rank_score","mean_percent_identity","mean_aligned_coverage_reference","mean_aligned_coverage_consensus","mean_orf_length_aa","reference_aa_length","orf_length_aa_range","total_internal_stop_count","max_internal_stop_count","mean_ambiguous_fraction","mean_covered_fraction","mean_depth","selection_reason","consensus_fastas"]
    if best_rank_tsv.exists():
        src = pd.read_csv(best_rank_tsv, sep="\t")
    else:
        src = targeted_df[targeted_df["target_type"] == "missing_gene_candidate"].copy()
        if src.empty:
            pd.DataFrame(columns=cols).to_csv(cross_tsv, sep="\t", index=False)
            cross_md.write_text("# Cross-readset missing-gene candidates\n\nNo missing-gene candidates.\n", encoding="utf-8")
            return cross_tsv, cross_md
        src["rank"] = 1
        src["rank_score"] = 0
    for c in ["rank","rank_score","percent_identity","aligned_coverage_reference","aligned_coverage_consensus","orf_length_aa","reference_aa_length","internal_stop_count","ambiguous_fraction","covered_fraction","mean_depth"]:
        if c in src.columns:
            src[c] = pd.to_numeric(src[c], errors="coerce")
    rows = []
    grp = src.groupby(["gene","target_id"], dropna=False)
    for (g, tid), d in grp:
        rs = sorted(set(d["read_set"].astype(str)))
        nrs = len(rs)
        best_rank = d["rank"].min() if "rank" in d.columns else 1
        mean_rank = d["rank"].mean() if "rank" in d.columns else 1
        mean_score = d["rank_score"].mean() if "rank_score" in d.columns else 0
        mean_pid = d["percent_identity"].mean()
        mean_covr = d["aligned_coverage_reference"].mean()
        mean_covc = d["aligned_coverage_consensus"].mean()
        mean_orf = d["orf_length_aa"].mean()
        ref_aa = d["reference_aa_length"].dropna().iloc[0] if d["reference_aa_length"].notna().any() else None
        min_orf = d["orf_length_aa"].min(); max_orf = d["orf_length_aa"].max()
        tot_stops = d["internal_stop_count"].fillna(0).sum()
        max_stops = d["internal_stop_count"].fillna(0).max()
        mean_amb = d["ambiguous_fraction"].mean()
        mean_covf = d["covered_fraction"].mean()
        mean_depth = d["mean_depth"].mean()
        if nrs >= 2 and tot_stops == 0 and mean_covr >= 50 and mean_pid >= 30:
            rec = "CROSS_READSET_GENE_CANDIDATE_SUPPORTED"
        elif nrs >= 2 and tot_stops == 0 and mean_covr >= 25 and mean_pid >= 25:
            rec = "CROSS_READSET_PARTIAL_GENE_CANDIDATE_SUPPORTED"
        elif nrs == 1 and max_stops == 0 and mean_covr >= 50 and mean_pid >= 30:
            rec = "SINGLE_READSET_GENE_CANDIDATE_SUPPORTED"
        else:
            rec = "MANUAL_REVIEW"
        consensus_fastas = ";".join(sorted(set(d.get("consensus_fasta", []).astype(str))))
        rows.append({"gene": g, "target_id": tid, "combined_recommendation": rec, "read_sets_supporting": ",".join(rs), "n_read_sets": nrs, "best_rank": best_rank, "mean_rank": round(mean_rank,3), "mean_rank_score": round(mean_score,3), "mean_percent_identity": round(mean_pid,3), "mean_aligned_coverage_reference": round(mean_covr,3), "mean_aligned_coverage_consensus": round(mean_covc,3), "mean_orf_length_aa": round(mean_orf,3), "reference_aa_length": ref_aa if ref_aa is not None else ".", "orf_length_aa_range": f"{min_orf}-{max_orf}", "total_internal_stop_count": int(tot_stops), "max_internal_stop_count": int(max_stops), "mean_ambiguous_fraction": round(mean_amb,4), "mean_covered_fraction": round(mean_covf,4), "mean_depth": round(mean_depth,3), "selection_reason": "n_read_sets+rank+coverage+identity+length+ambiguity+depth", "consensus_fastas": consensus_fastas})
    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(columns=cols)
    else:
        out = out.sort_values(["n_read_sets","best_rank","mean_rank","mean_aligned_coverage_reference","mean_percent_identity","mean_ambiguous_fraction","mean_covered_fraction","mean_depth"], ascending=[False,True,True,False,False,True,False,False]).reset_index(drop=True)
        out["combined_rank"] = range(1, len(out) + 1)
        out = out[cols]
    out.to_csv(cross_tsv, sep="\t", index=False)
    with open(cross_md, "w", encoding="utf-8") as md:
        md.write("# Cross-readset missing-gene candidates\n\n")
        if out.empty:
            md.write("No missing-gene candidates.\n")
        else:
            for gene, gdf in out.groupby("gene"):
                md.write(f"## {gene}\n\n")
                best = gdf.iloc[0]
                md.write("Best candidate:\n")
                md.write(f"- target_id: {best['target_id']}\n")
                md.write(f"- combined_recommendation: {best['combined_recommendation']}\n")
                md.write(f"- read_sets_supporting: {best['read_sets_supporting']}\n")
                md.write(f"- mean identity: {best['mean_percent_identity']}\n")
                md.write(f"- mean reference coverage: {best['mean_aligned_coverage_reference']}\n")
                md.write(f"- ORF length range: {best['orf_length_aa_range']}\n\n")
                md.write("Ranked candidates:\n")
                for rr in gdf.itertuples():
                    md.write(f"- rank {rr.combined_rank}: {rr.target_id} ({rr.combined_recommendation}) [{rr.read_sets_supporting}]\n")
                md.write("\n")
        md.write("> Nenhum GenBank foi alterado.\n")
    return cross_tsv, cross_md
