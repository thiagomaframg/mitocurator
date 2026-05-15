from __future__ import annotations
from pathlib import Path
import gzip
import random
import shutil
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

from .utils import ensure_dir, safe_get, run_cmd, get_genetic_code
from .io import read_record
from .targeted_consensus import _reference_protein

LONG_TYPES = {"hifi", "pacbio_hifi", "pacbio_clr", "clr", "ont", "nanopore"}
SHORT_TYPES = {"illumina", "illumina_pe", "pe", "illumina_se", "se"}


def _type_from_readset(config, read_set):
    for grp in (safe_get(config, ["reads", "long"], []) or []) + (safe_get(config, ["reads", "short"], []) or []):
        if isinstance(grp, dict) and str(grp.get("name")) == read_set:
            return str(grp.get("type", ""))
    for rs in (safe_get(config, ["read_support", "read_sets"], []) or []):
        if str(rs.get("name")) == read_set:
            return str(rs.get("type", ""))
    n = read_set.lower()
    if "illumina" in n:
        return "illumina_pe"
    if "hifi" in n:
        return "pacbio_hifi"
    if "ont" in n:
        return "ont"
    return ""


def _read_base_id(name: str) -> str:
    nid = (name or "").split()[0]
    if nid.endswith("/1") or nid.endswith("/2"):
        nid = nid[:-2]
    return nid


def _pick_ids(ids, limit, rng):
    if limit is None or limit < 0 or len(ids) <= limit:
        return sorted(ids)
    return sorted(rng.sample(list(ids), limit))


def _collect_single_records(fastq_path: Path):
    recs = {}
    if not fastq_path or str(fastq_path) == "." or not Path(fastq_path).exists():
        return recs
    with gzip.open(fastq_path, "rt") if str(fastq_path).endswith(".gz") else open(fastq_path, "rt", encoding="utf-8") as h:
        for rec in SeqIO.parse(h, "fastq"):
            recs[_read_base_id(rec.id)] = rec
    return recs


def _collect_paired_records(r1_path: Path, r2_path: Path):
    pairs = {}
    if any((not p, str(p) == ".", not Path(p).exists()) for p in [r1_path, r2_path]):
        return pairs
    op1 = gzip.open(r1_path, "rt") if str(r1_path).endswith(".gz") else open(r1_path, "rt", encoding="utf-8")
    op2 = gzip.open(r2_path, "rt") if str(r2_path).endswith(".gz") else open(r2_path, "rt", encoding="utf-8")
    with op1 as h1, op2 as h2:
        for a, b in zip(SeqIO.parse(h1, "fastq"), SeqIO.parse(h2, "fastq")):
            ia, ib = _read_base_id(a.id), _read_base_id(b.id)
            if ia == ib:
                pairs[ia] = (a, b)
    return pairs




def _scan_orfs(seq: str, code: int, min_orf_nt: int = 150):
    cands = []
    strands = [("+", seq), ("-", str(Seq(seq).reverse_complement()))]
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
                    cands.append({"orf_id": f"orf{oid}", "aa": aa, "internal_stop_count": aa[:-1].count("*")})
                    oid += 1
                i = j + 3
    return cands

def _align_metrics(query_aa: str, ref_aa: str):
    if not query_aa or not ref_aa:
        return ".", ".", "."
    al = PairwiseAligner(mode="global")
    aln = al.align(query_aa, ref_aa)[0]
    q, r = aln[0], aln[1]
    matches = sum(1 for a, b in zip(q, r) if a == b and a != "-" and b != "-")
    aligned = sum(1 for a, b in zip(q, r) if a != "-" and b != "-")
    pid = round(100.0 * matches / max(aligned, 1), 2)
    covr = round(100.0 * aligned / max(len(ref_aa), 1), 2)
    return pid, covr, float(aln.score)


def _load_target_regions(targeted_dir: Path):
    """Load target regions robustly from targets.bed or targeted_read_extraction.tsv.

    Supports:
      1) BED-like: contig start0 end target_id
      2) target-first: target_id contig start0 end

    Skips headers such as: target_id, contig, start0, end.
    """
    regions = {}

    def _is_header(parts):
        low = [str(x).strip().lower() for x in parts]
        header_tokens = {
            "target_id", "target", "gene",
            "contig", "chrom", "seqid", "sequence",
            "start", "start0", "start_0",
            "end", "stop"
        }
        return any(x in header_tokens for x in low[:4])

    def _add_region(parts):
        if len(parts) < 4:
            return
        if _is_header(parts):
            return

        # Format A: contig start0 end target_id
        try:
            contig = str(parts[0]).strip()
            start0 = int(parts[1])
            end = int(parts[2])
            target_id = str(parts[3]).strip()
            regions[target_id] = (contig, start0, end)
            return
        except Exception:
            pass

        # Format B: target_id contig start0 end
        try:
            target_id = str(parts[0]).strip()
            contig = str(parts[1]).strip()
            start0 = int(parts[2])
            end = int(parts[3])
            regions[target_id] = (contig, start0, end)
            return
        except Exception:
            return

    bed = targeted_dir / "targets.bed"
    if bed.exists():
        with bed.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                _add_region(line.split("\t"))

    # Fallback: targeted_read_extraction.tsv
    tsv = targeted_dir / "targeted_read_extraction.tsv"
    if tsv.exists():
        try:
            import pandas as pd
            df = pd.read_csv(tsv, sep="\t")
            cols = {c.lower(): c for c in df.columns}

            tid_col = cols.get("target_id") or cols.get("target")
            contig_col = cols.get("contig") or cols.get("seqid") or cols.get("chrom")
            start_col = cols.get("start0") or cols.get("start") or cols.get("start_0")
            end_col = cols.get("end") or cols.get("stop")

            if tid_col and contig_col and start_col and end_col:
                for _, row in df.iterrows():
                    tid = str(row[tid_col])
                    if tid in regions:
                        continue
                    try:
                        regions[tid] = (
                            str(row[contig_col]),
                            int(row[start_col]),
                            int(row[end_col]),
                        )
                    except Exception:
                        continue
        except Exception:
            pass

    return regions


def _select_read_ids_from_bam(bam_path: Path, region, max_reads: int, min_mapping_quality: int, seed_text: str, source_type: str):
    rows = []
    if not bam_path.exists():
        return rows
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        it = bam.fetch(region[0], max(0, int(region[1])), int(region[2])) if region else bam.fetch(until_eof=True)
        for r in it:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if int(r.mapping_quality) < int(min_mapping_quality):
                continue
            rid = _read_base_id(r.query_name)
            rows.append({"read_id": rid, "source_type": source_type, "mapping_contig": bam.get_reference_name(r.reference_id) if r.reference_id >= 0 else ".", "mapping_start": int(r.reference_start)+1, "mapping_end": int(r.reference_end or r.reference_start), "mapping_quality": int(r.mapping_quality), "read_length": int(r.query_length or 0)})
    rows = sorted(rows, key=lambda x: x["read_id"])
    uniq = {}
    for r in rows:
        uniq[r["read_id"]] = r
    rows = [uniq[k] for k in sorted(uniq.keys())]
    rng = random.Random(seed_text)
    if max_reads >= 0 and len(rows) > max_reads:
        idx = sorted(rng.sample(range(len(rows)), max_reads))
        rows = [rows[i] for i in idx]
    for i, r in enumerate(rows, 1):
        r["selected_order"] = i
    return rows


def normalize_read_id(header_or_id: str) -> str:
    x = str(header_or_id or "").strip()
    if x.startswith("@"):
        x = x[1:]
    x = x.split()[0]
    if x.endswith("/1") or x.endswith("/2"):
        x = x[:-2]
    return x


def _iter_fastq(path: Path):
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8") as h:
        yield from SeqIO.parse(h, "fastq")


def _as_list(v):
    if v is None: return []
    if isinstance(v, list): return [str(x) for x in v]
    return [str(v)]


def _recover_single_from_sources(selected_ids, paths, out_fastq):
    found = {}
    for p in [Path(x) for x in paths if x and x != "."]:
        if not p.exists():
            continue
        for rec in _iter_fastq(p):
            rid = normalize_read_id(rec.description or rec.id)
            if rid in selected_ids and rid not in found:
                found[rid] = rec
                if len(found) == len(selected_ids):
                    break
    ordered = [rid for rid in selected_ids if rid in found]
    if ordered:
        with gzip.open(out_fastq, "wt") as oh:
            SeqIO.write([found[r] for r in ordered], oh, "fastq")
    return len(ordered), len(selected_ids) - len(ordered), sum(len(found[r].seq) for r in ordered)

def _recover_pe(selected_ids, r1_path, r2_path, out_r1, out_r2):
    found = {}
    if not (Path(r1_path).exists() and Path(r2_path).exists()):
        return 0, len(selected_ids), 0
    for a,b in zip(_iter_fastq(Path(r1_path)), _iter_fastq(Path(r2_path))):
        i1 = normalize_read_id(a.description or a.id); i2 = normalize_read_id(b.description or b.id)
        if i1 != i2:
            continue
        if i1 in selected_ids and i1 not in found:
            found[i1]=(a,b)
            if len(found)==len(selected_ids):
                break
    ordered=[rid for rid in selected_ids if rid in found]
    if ordered:
        with gzip.open(out_r1,'wt') as h1, gzip.open(out_r2,'wt') as h2:
            SeqIO.write([found[r][0] for r in ordered], h1, 'fastq')
            SeqIO.write([found[r][1] for r in ordered], h2, 'fastq')
    return len(ordered), len(selected_ids)-len(ordered), sum(len(found[r][0].seq)+len(found[r][1].seq) for r in ordered)


def _write_fastq_from_bam_selected(bam_path: Path, selected_ids, out_fastq: Path, paired=False, out_r1: Path=None, out_r2: Path=None):
    selected = set(selected_ids)
    if not bam_path.exists() or not selected:
        return 0, 0, 0
    seen = set()
    total_bases = 0
    lengths = []
    if paired:
        pairs = {}
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for r in bam.fetch(until_eof=True):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = normalize_read_id(r.query_name)
                if rid not in selected or not r.query_sequence:
                    continue
                q = r.query_qualities or [40] * len(r.query_sequence)
                qual = ''.join(chr(min(max(int(x),0),93)+33) for x in q)
                ent = (r.query_name, r.query_sequence, qual)
                d = pairs.setdefault(rid, {})
                d['r1' if r.is_read1 else 'r2' if r.is_read2 else 'r1'] = ent
        written=0
        with gzip.open(out_r1,'wt') as h1, gzip.open(out_r2,'wt') as h2:
            for rid in selected_ids:
                d=pairs.get(rid,{})
                if 'r1' in d and 'r2' in d:
                    h,seq,qual = d['r1']; h1.write("@%s\n%s\n+\n%s\n" % (h,seq,qual))
                    h,seq,qual = d['r2']; h2.write("@%s\n%s\n+\n%s\n" % (h,seq,qual))
                    written += 1
                    total_bases += len(d['r1'][1]) + len(d['r2'][1]); lengths.append(len(d['r1'][1])+len(d['r2'][1]))
        return written, len(selected_ids)-written, total_bases
    else:
        recs = {}
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for r in bam.fetch(until_eof=True):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = normalize_read_id(r.query_name)
                if rid in selected and rid not in recs and r.query_sequence:
                    q = r.query_qualities or [40] * len(r.query_sequence)
                    qual = ''.join(chr(min(max(int(x),0),93)+33) for x in q)
                    recs[rid] = (r.query_name, r.query_sequence, qual)
        with gzip.open(out_fastq,'wt') as oh:
            for rid in selected_ids:
                if rid in recs:
                    h,seq,qual = recs[rid]
                    oh.write("@%s\n%s\n+\n%s\n" % (h,seq,qual))
                    total_bases += len(seq); lengths.append(len(seq))
        return len([x for x in selected_ids if x in recs]), len([x for x in selected_ids if x not in recs]), total_bases

def _merged_coverage(intervals, qlen):
    if not intervals or qlen <= 0:
        return 0.0
    iv=sorted((min(a,b),max(a,b)) for a,b in intervals)
    merged=[]
    for a,b in iv:
        if not merged or a>merged[-1][1]+1: merged.append([a,b])
        else: merged[-1][1]=max(merged[-1][1],b)
    cov=sum(b-a+1 for a,b in merged)
    return round(min(100.0,100.0*cov/max(qlen,1)),2)

def run_candidate_assembly(config, root: Path, consensus_dir: Path, pools_dir: Path, refinement_dir: Path, outdir: Path):
    outdir = ensure_dir(outdir)
    targets_tsv = outdir / "candidate_assembly_targets.tsv"
    sum_tsv = outdir / "candidate_assembly_summary.tsv"
    sum_md = outdir / "candidate_assembly_summary.md"
    ca = safe_get(config, ["candidate_assembly"], {}) or {}
    strategy = str(ca.get("assembly_pool_strategy", "mapped_mitogenome_plus_target"))
    target_req = int(ca.get("max_target_reads_per_candidate", 85))
    mito_req = int(ca.get("max_mitogenome_reads_per_candidate", 300))
    rng = random.Random(int(ca.get("random_seed", 42)))

    src_mode = ca.get("candidate_source", "cross_readset")
    cross = consensus_dir / "cross_readset_missing_gene_candidates.tsv"
    best = consensus_dir / "best_missing_gene_candidates.tsv"
    if src_mode == "cross_readset" and cross.exists():
        cand = pd.read_csv(cross, sep="\t").rename(columns={"combined_rank": "rank", "combined_recommendation": "recommendation"})
    elif best.exists():
        cand = pd.read_csv(best, sep="\t").rename(columns={"rank": "rank", "recommendation": "recommendation"})
    else:
        pd.DataFrame().to_csv(targets_tsv, sep="\t", index=False)
        pd.DataFrame().to_csv(sum_tsv, sep="\t", index=False)
        sum_md.write_text("# Candidate assembly\n\nNo candidate source found.\n", encoding="utf-8")
        return outdir

    inc = set(ca.get("include_recommendations", ["CROSS_READSET_GENE_CANDIDATE_SUPPORTED", "CROSS_READSET_PARTIAL_GENE_CANDIDATE_SUPPORTED", "SINGLE_READSET_GENE_CANDIDATE_SUPPORTED"]))
    max_per_gene = int(ca.get("max_candidates_per_gene", 1))
    cand = cand[cand["recommendation"].isin(inc)]
    if "rank" in cand.columns:
        cand = cand[cand["rank"] <= max_per_gene]

    pools = pd.read_csv(pools_dir / "reconstruction_pools.tsv", sep="\t") if (pools_dir / "reconstruction_pools.tsv").exists() else pd.DataFrame()
    rec, _ = read_record(refinement_dir / "refined.gb")
    refined_fa = outdir / "refined_mitogenome.fasta"
    SeqIO.write([SeqIO.SeqRecord(rec.seq, id=rec.id, description="")], str(refined_fa), "fasta")
    run_cmd(["makeblastdb", "-in", str(refined_fa), "-dbtype", "nucl"], check=False)
    threads = int(ca.get("threads", 16))
    min_len = int(ca.get("min_contig_len", 300))
    code = get_genetic_code(config, default=5)
    ps = safe_get(config, ["candidate_assembly", "protein_search"], {}) or {}
    p_word = int(ps.get("word_size", 2))
    p_matrix = str(ps.get("matrix", "PAM30"))
    p_comp = int(ps.get("comp_based_stats", 0))
    p_seg = "no" if bool(ps.get("disable_seg", True)) else "yes"
    p_soft = "false" if not bool(ps.get("soft_masking", False)) else "true"
    run_blastx_fallback = bool(ps.get("run_blastx_fallback", True))
    run_orfscan_fallback = bool(ps.get("run_orfscan_fallback", True))
    min_orf_nt = int(ps.get("min_orf_nt", 150))

    rows = []
    target_rows = []
    target_regions = _load_target_regions(root / "11_targeted_extraction")
    for rr in cand.itertuples():
        gene = str(rr.gene)
        tid = str(rr.target_id)
        base = ensure_dir(outdir / gene / tid)
        ad = ensure_dir(base / "assembly")
        bd = ensure_dir(base / "blast")
        dd = ensure_dir(base / "diagnosis")
        rd = ensure_dir(base / "reads_downsampled")
        subset = pools[(pools["target_id"].astype(str) == tid)]
        read_sets = sorted(set(subset["read_set"].astype(str)))

        for rs in read_sets:
            rs_sub = subset[subset["read_set"].astype(str) == rs]
            target_only = rs_sub[rs_sub["pool_type"].astype(str) == "target_only"]
            mito_only = rs_sub[rs_sub["pool_type"].astype(str) == "mitogenome_mapped"]
            combined = rs_sub[rs_sub["pool_type"].astype(str) == "combined"]
            rtype = _type_from_readset(config, rs)
            assembler = None
            asm_fa = None
            status = "ASSEMBLY_NOT_RUN"
            in_fastq = in_r1 = in_r2 = "."

            t_avail = m_avail = t_used = m_used = dedup_rm = final_used = 0
            selected_ids_tsv = "."
            total_bases = 0
            n50 = 0
            if strategy == "mapped_mitogenome_plus_target":
                bam = root / "08_read_mapping" / f"{rs}.sorted.bam"
                region = target_regions.get(tid)
                seed_base = f"{ca.get('random_seed',42)}|{gene}|{tid}|{rs}"
                mrows = _select_read_ids_from_bam(bam, None, mito_req, int(ca.get('min_mapping_quality',20)), seed_base+"|mitogenome_mapped", "mitogenome_mapped")
                trows = _select_read_ids_from_bam(bam, region, target_req, int(ca.get('min_mapping_quality',20)), seed_base+"|target_mapped", "target_mapped")
                t_avail, m_avail = len(trows), len(mrows)
                merged = {}
                for rrw in trows + mrows:
                    if rrw["read_id"] not in merged:
                        merged[rrw["read_id"]] = rrw
                dedup_rm = len(trows) + len(mrows) - len(merged)
                t_used = len([x for x in merged.values() if x["source_type"] == "target_mapped"])
                m_used = len([x for x in merged.values() if x["source_type"] == "mitogenome_mapped"])
                final_used = len(merged)
                if merged:
                    selected_ids_tsv = str(rd / f"{tid}.{rs}.selected_read_ids.tsv")
                    pd.DataFrame(sorted(merged.values(), key=lambda x: x["selected_order"] if "selected_order" in x else 0)).to_csv(selected_ids_tsv, sep="\t", index=False)
                sel_ids = [x["read_id"] for x in sorted(merged.values(), key=lambda x: x.get("selected_order", 0))]
                if sel_ids:
                    use_bam_fastq = bool(ca.get("use_bam_sequences_for_fastq", True))
                    if rtype in LONG_TYPES or rtype in {"illumina_se", "se"}:
                        in_fastq = str(rd / f"{tid}.{rs}.mito300_target85.fastq.gz")
                        recov, miss, total_bases = _write_fastq_from_bam_selected(bam, sel_ids, Path(in_fastq), paired=False)
                        final_used = recov
                        print(f"[candidate_assembly] read_set={rs} selected={len(sel_ids)} recovered={recov} missing={miss} fastq={in_fastq}")
                    else:
                        in_r1 = str(rd / f"{tid}.{rs}.mito_cov250_target_cov300_R1.fastq.gz")
                        in_r2 = str(rd / f"{tid}.{rs}.mito_cov250_target_cov300_R2.fastq.gz")
                        recov, miss, total_bases = _write_fastq_from_bam_selected(bam, sel_ids, None, paired=True, out_r1=Path(in_r1), out_r2=Path(in_r2))
                        final_used = recov
                        print(f"[candidate_assembly] read_set={rs} selected_pairs={len(sel_ids)} recovered_pairs={recov} missing_pairs={miss} R1={in_r1} R2={in_r2}")
                    if final_used > 0:
                        lens = [int(x.get("read_length",0)) for x in merged.values() if int(x.get("read_length",0))>0]
                        if lens:
                            a=sorted(lens, reverse=True); c=0; half=sum(a)/2
                            for L in a:
                                c += L
                                if c >= half:
                                    n50 = L; break
                    elif len(sel_ids) > 0:
                        status = "SELECTED_READS_NOT_FOUND_IN_FASTQ"

            elif strategy == "targeted_plus_mitogenome":
                if rtype in LONG_TYPES or rtype in {"illumina_se", "se"}:
                    tref = target_only.iloc[0] if not target_only.empty else None
                    mref = mito_only.iloc[0] if not mito_only.empty else None
                    t_records = _collect_single_records(Path(str(tref.output_fastq))) if tref is not None else {}
                    m_records = _collect_single_records(Path(str(mref.output_fastq))) if mref is not None else {}
                    t_ids, m_ids = set(t_records.keys()), set(m_records.keys())
                    t_avail, m_avail = len(t_ids), len(m_ids)
                    sel_t = _pick_ids(t_ids, target_req, rng)
                    sel_m = _pick_ids(m_ids, mito_req, rng)
                    chosen = {}
                    for i in sel_t:
                        chosen[i] = t_records[i]
                    for i in sel_m:
                        if i not in chosen:
                            chosen[i] = m_records[i]
                    t_used = len(set(sel_t))
                    m_used = len([i for i in sel_m if i not in set(sel_t)])
                    dedup_rm = len(sel_t) + len(sel_m) - len(chosen)
                    final_used = len(chosen)
                    if chosen:
                        in_fastq = str(rd / f"{tid}.{rs}.target85_mito300.fastq.gz")
                        with gzip.open(in_fastq, "wt") as oh:
                            SeqIO.write([chosen[i] for i in sorted(chosen.keys())], oh, "fastq")
                elif rtype in SHORT_TYPES:
                    pref_t = target_only.iloc[0] if not target_only.empty else None
                    pref_m = mito_only.iloc[0] if not mito_only.empty else None
                    t_pairs = _collect_paired_records(Path(str(pref_t.output_fastq_r1)), Path(str(pref_t.output_fastq_r2))) if pref_t is not None and str(pref_t.output_format) == "paired_fastq" else {}
                    m_pairs = _collect_paired_records(Path(str(pref_m.output_fastq_r1)), Path(str(pref_m.output_fastq_r2))) if pref_m is not None and str(pref_m.output_format) == "paired_fastq" else {}
                    t_ids, m_ids = set(t_pairs.keys()), set(m_pairs.keys())
                    t_avail, m_avail = len(t_ids), len(m_ids)
                    sel_t = _pick_ids(t_ids, target_req, rng)
                    sel_m = _pick_ids(m_ids, mito_req, rng)
                    chosen = {}
                    for i in sel_t:
                        chosen[i] = t_pairs[i]
                    for i in sel_m:
                        if i not in chosen:
                            chosen[i] = m_pairs[i]
                    t_used = len(set(sel_t))
                    m_used = len([i for i in sel_m if i not in set(sel_t)])
                    dedup_rm = len(sel_t) + len(sel_m) - len(chosen)
                    final_used = len(chosen)
                    if chosen:
                        in_r1 = str(rd / f"{tid}.{rs}.target85_mito300_R1.fastq.gz")
                        in_r2 = str(rd / f"{tid}.{rs}.target85_mito300_R2.fastq.gz")
                        with gzip.open(in_r1, "wt") as h1, gzip.open(in_r2, "wt") as h2:
                            SeqIO.write([chosen[i][0] for i in sorted(chosen.keys())], h1, "fastq")
                            SeqIO.write([chosen[i][1] for i in sorted(chosen.keys())], h2, "fastq")

            if rtype in LONG_TYPES and bool(ca.get("run_long_read_assembly", True)):
                assembler = str(ca.get("long_read_assembler", "flye"))
                if in_fastq != ".":
                    mode = "--pacbio-hifi" if "hifi" in rtype else "--pacbio-raw" if "clr" in rtype else "--nano-hq"
                    outd = ad / f"flye_{rs}"
                    run_cmd(["bash", "-lc", f"flye {mode} {in_fastq} --out-dir {outd} --threads {threads} --genome-size {safe_get(config,['candidate_assembly','flye','genome_size'],'20k')} {safe_get(config,['candidate_assembly','flye','extra_args'],'')}"], check=False)
                    asm = outd / "assembly.fasta"
                    if asm.exists():
                        asm_fa = ad / f"{tid}.{rs}.flye.fasta"
                        shutil.copyfile(asm, asm_fa)
                        status = "OK"
                    else:
                        status = "ASSEMBLY_FAILED"
                else:
                    status = "NO_INPUT_READS"
            elif rtype in SHORT_TYPES and bool(ca.get("run_short_read_assembly", True)):
                assembler = str(ca.get("short_read_assembler", "spades"))
                outd = ad / f"spades_{rs}"
                if in_r1 != "." and in_r2 != ".":
                    cmd = f"spades.py -1 {in_r1} -2 {in_r2} -o {outd} --threads {threads} {safe_get(config,['candidate_assembly','spades','extra_args'],'')}"
                elif in_fastq != ".":
                    cmd = f"spades.py -s {in_fastq} -o {outd} --threads {threads} {safe_get(config,['candidate_assembly','spades','extra_args'],'')}"
                else:
                    cmd = None
                if cmd:
                    run_cmd(["bash", "-lc", cmd], check=False)
                    asm = outd / "contigs.fasta"
                    if asm.exists():
                        asm_fa = ad / f"{tid}.{rs}.spades.fasta"
                        shutil.copyfile(asm, asm_fa)
                        status = "OK"
                    else:
                        status = "ASSEMBLY_FAILED"
                else:
                    status = "NO_INPUT_READS"
            else:
                status = "ASSEMBLER_NOT_AVAILABLE"

            best_contig = "."
            best_len = 0
            bident = bqcov = bsstart = bsend = "."
            tbpid = tbqcov = tbscov = tbbs = "."
            recm = "NO_LOCAL_ASSEMBLY_SUPPORT"
            protein_method = "."; tblastn_status = blastx_status = orfscan_status = "NOT_RUN"
            protein_best_pid = protein_best_cov = protein_best_score = "."
            protein_num_hits = 0; protein_best_hsp_pident = "."; protein_best_hsp_bitscore = "."; protein_best_contig_id = "."
            if asm_fa and Path(asm_fa).exists():
                sel = ad / "selected_candidate_contigs.fasta"
                kept = [s for s in SeqIO.parse(str(asm_fa), "fasta") if len(s.seq) >= min_len]
                if kept:
                    SeqIO.write(kept, str(sel), "fasta")
                    k = max(kept, key=lambda x: len(x.seq))
                    best_len, best_contig = len(k.seq), k.id
                    b6 = bd / "blastn_vs_refined_mitogenome.tsv"
                    run_cmd(["blastn", "-query", str(sel), "-db", str(refined_fa), "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen", "-evalue", str(ca.get("blastn_evalue", "1e-10")), "-num_threads", str(threads), "-out", str(b6)], check=False)
                    if b6.exists() and b6.stat().st_size > 0:
                        x = pd.read_csv(b6, sep="\t", header=None).iloc[0]
                        bident = float(x[2]); qlen = float(x[12]); bqcov = round(100 * float(x[3]) / max(qlen, 1), 2); bsstart = int(x[8]); bsend = int(x[9])
                    refaa = _reference_protein(config, gene, code)
                    if refaa:
                        qfaa = bd / f"reference_{gene}.faa"
                        with open(qfaa, "w", encoding="utf-8") as f:
                            f.write(f">{gene}\n{refaa}\n")
                        run_cmd(["makeblastdb", "-in", str(sel), "-dbtype", "nucl"], check=False)
                        t6 = bd / "tblastn_reference_protein_vs_candidate_assembly.tsv"
                        run_cmd(["tblastn", "-query", str(qfaa), "-db", str(sel), "-db_gencode", str(code), "-seg", p_seg, "-soft_masking", p_soft, "-word_size", str(p_word), "-matrix", p_matrix, "-comp_based_stats", str(p_comp), "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen", "-evalue", str(ca.get("tblastn_evalue", "1e-5")), "-num_threads", str(threads), "-out", str(t6)], check=False)
                        if t6.exists() and t6.stat().st_size > 0:
                            y = pd.read_csv(t6, sep="\t", header=None).iloc[0]
                            tbpid = float(y[2]); qlen = float(y[12]); slen = float(y[13]); al = float(y[3]); tbqcov = round(100 * al / max(qlen, 1), 2); tbscov = round(100 * al / max(slen, 1), 2); tbbs = float(y[11])
                            protein_method = "tblastn"; tblastn_status = "HIT"; protein_best_pid = tbpid; protein_best_cov = tbqcov; protein_best_score = tbbs
                        else:
                            tblastn_status = "NO_HIT"
                            if run_blastx_fallback:
                                pdb = bd / f"reference_{gene}_prot_db.faa"
                                with open(pdb, "w", encoding="utf-8") as f:
                                    f.write(f">{gene}\n{refaa}\n")
                                run_cmd(["makeblastdb", "-in", str(pdb), "-dbtype", "prot"], check=False)
                                x6 = bd / "blastx_reference_protein_vs_candidate_assembly.tsv"
                                run_cmd(["blastx", "-query", str(sel), "-db", str(pdb), "-query_gencode", str(code), "-seg", p_seg, "-soft_masking", p_soft, "-word_size", str(p_word), "-matrix", p_matrix, "-comp_based_stats", str(p_comp), "-evalue", str(ca.get("tblastn_evalue", "1e-5")), "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen", "-num_threads", str(threads), "-out", str(x6)], check=False)
                                if x6.exists() and x6.stat().st_size > 0:
                                    y = pd.read_csv(x6, sep="\t", header=None).iloc[0]
                                    protein_method = "blastx"; blastx_status = "HIT"; protein_best_pid = float(y[2]); protein_best_cov = round(100*float(y[3])/max(float(y[12]),1),2); protein_best_score = float(y[11])
                                else:
                                    blastx_status = "NO_HIT"
                            if protein_method == "." and run_orfscan_fallback:
                                best_orf = None
                                for contig in kept:
                                    for orf in _scan_orfs(str(contig.seq).upper(), code, min_orf_nt=min_orf_nt):
                                        pid, covr, score = _align_metrics(orf["aa"], refaa)
                                        if pid == ".":
                                            continue
                                        cand = (float(score), float(pid), float(covr), -int(orf["internal_stop_count"]))
                                        if best_orf is None or cand > best_orf[0]:
                                            best_orf = (cand, pid, covr, score, orf["internal_stop_count"])
                                if best_orf:
                                    protein_method = "orfscan_pairwise"; orfscan_status = "HIT"
                                    protein_best_pid, protein_best_cov, protein_best_score = best_orf[1], best_orf[2], best_orf[3]
                                    pd.DataFrame([{"method":"orfscan_pairwise","pident":protein_best_pid,"coverage_reference":protein_best_cov,"score":protein_best_score,"internal_stop_count":best_orf[4]}]).to_csv(bd / "orfscan_reference_protein_vs_candidate_assembly.tsv", sep="\t", index=False)
                                else:
                                    orfscan_status = "NO_HIT"
                    if protein_method != "." and protein_best_pid != "." and float(protein_best_pid) >= float(ca.get("tblastn_min_identity", 25)) and float(protein_best_cov) >= float(ca.get("tblastn_min_coverage", 25)):
                        recm = "LOCAL_ASSEMBLY_SUPPORTS_MISSING_GENE" if protein_method != "orfscan_pairwise" else "LOCAL_ASSEMBLY_SUPPORTS_MISSING_GENE_BY_ORF_ALIGNMENT"
                        if bqcov != "." and bqcov < float(ca.get("blastn_min_coverage", 30)):
                            recm = "LOCAL_ASSEMBLY_SUGGESTS_MITOGENOME_ERROR"
                    elif bident != ".":
                        recm = "LOCAL_ASSEMBLY_MITOCHONDRIAL_BUT_NO_PROTEIN_HIT"
            row = {
                "gene": gene, "target_id": tid, "read_set": rs,
                "assembly_pool_strategy": strategy,
                "target_reads_requested": target_req, "target_reads_available": t_avail, "target_reads_used": t_used,
                "mitogenome_reads_requested": mito_req, "mitogenome_reads_available": m_avail, "mitogenome_reads_used": m_used,
                "duplicated_reads_removed": dedup_rm, "final_reads_used_for_assembly": final_used,
                "reads_recovered_from_fastq": final_used if in_fastq != "." else 0, "read_pairs_recovered_from_fastq": final_used if in_r1 != "." else 0, "reads_missing_in_fastq": max(0,(t_used+m_used)-final_used), "fastq_source": "BAM_QUERY_SEQUENCE", "selected_read_ids_tsv": selected_ids_tsv, "total_selected_bases": total_bases, "selected_read_n50": n50, "assembly_input_fastq": in_fastq, "assembly_input_r1": in_r1, "assembly_input_r2": in_r2,
                "assembler": assembler or ".", "assembly_fasta": str(asm_fa) if asm_fa else ".",
                "best_contig_id": best_contig, "best_contig_len": best_len,
                "blastn_best_hit_refined_seqid": rec.id if bident != "." else ".", "blastn_best_pident": bident,
                "blastn_best_qcov": bqcov, "blastn_best_sstart": bsstart, "blastn_best_send": bsend,
                "tblastn_best_pident": tbpid, "tblastn_best_query_coverage": tbqcov, "tblastn_best_subject_coverage": tbscov,
                "tblastn_best_bitscore": tbbs, "reference_protein_length": len(_reference_protein(config, gene, code) or ""),
                "candidate_region_status": status, "recommendation": recm, "protein_search_method": protein_method, "tblastn_status": tblastn_status, "blastx_status": blastx_status, "orfscan_status": orfscan_status,
                "protein_search_best_pident": protein_best_pid, "protein_search_best_reference_coverage": protein_best_cov, "protein_search_best_bitscore_or_score": protein_best_score, "protein_search_num_hits": protein_num_hits, "protein_search_best_hsp_pident": protein_best_hsp_pident, "protein_search_best_hsp_bitscore": protein_best_hsp_bitscore, "protein_search_merged_reference_coverage": protein_best_cov, "protein_search_best_contig_id": protein_best_contig_id,
                "comment": "diagnostic-only",
            }
            rows.append(row)
            target_rows.append({k: row[k] for k in [
                "gene", "target_id", "read_set", "assembly_pool_strategy", "target_reads_requested", "target_reads_available", "target_reads_used",
                "mitogenome_reads_requested", "mitogenome_reads_available", "mitogenome_reads_used", "duplicated_reads_removed",
                "final_reads_used_for_assembly", "assembly_input_fastq", "assembly_input_r1", "assembly_input_r2"
            ]})
            pd.DataFrame([row]).to_csv(dd / "candidate_gene_diagnosis.tsv", sep="\t", index=False)

    pd.DataFrame(target_rows).to_csv(targets_tsv, sep="\t", index=False)
    sdf = pd.DataFrame(rows)
    sdf.to_csv(sum_tsv, sep="\t", index=False)
    with open(sum_md, "w", encoding="utf-8") as md:
        md.write("# Candidate assembly summary\n\n")
        if sdf.empty:
            md.write("No candidate assemblies generated.\n")
        else:
            for r in sdf.itertuples():
                md.write(f"## {r.gene} :: {r.target_id} :: {r.read_set}\n\n")
                md.write(f"- strategy: {r.assembly_pool_strategy}\n")
                md.write(f"- target_only usadas: {r.target_reads_used}/{r.target_reads_available} (req={r.target_reads_requested})\n")
                md.write(f"- mitogenome usadas: {r.mitogenome_reads_used}/{r.mitogenome_reads_available} (req={r.mitogenome_reads_requested})\n")
                md.write(f"- duplicatas removidas: {r.duplicated_reads_removed}\n")
                md.write(f"- reads finais para montagem: {r.final_reads_used_for_assembly}\n")
                md.write(f"- selected IDs TSV: {r.selected_read_ids_tsv}\n")
                md.write(f"- status montagem: {r.candidate_region_status}\n")
                md.write(f"- status BLASTN: {r.blastn_best_pident}\n")
                md.write(f"- status TBLASTN: {r.tblastn_status}\n")
                md.write(f"- status BLASTX: {r.blastx_status}\n")
                md.write(f"- status ORF scan: {r.orfscan_status}\n")
                md.write(f"- protein search: method={r.protein_search_method} tblastn={r.tblastn_status} blastx={r.blastx_status} orfscan={r.orfscan_status}\n")
                md.write(f"- protein best: pid={r.protein_search_best_pident} cov_ref={r.protein_search_best_reference_coverage} score={r.protein_search_best_bitscore_or_score}\n")
                md.write(f"- recomendação: {r.recommendation}\n\n")
        md.write("Nenhum GenBank foi alterado nesta etapa.\n")
    return outdir
