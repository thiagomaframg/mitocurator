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
    "illumina_pe": "sr", "illumina": "sr", "pe": "sr", "illumina_se": "sr",
}


def get_by_dotted_path(config: dict, path: str):
    cur = config
    for part in path.split("."):
        if isinstance(cur, list):
            try:
                cur = cur[int(part)]
                continue
            except Exception:
                return None
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def _normalize_reads_for_type(reads, rtype: str):
    if rtype in {"illumina_pe", "illumina", "pe"}:
        if isinstance(reads, dict):
            r1 = reads.get("r1") or reads.get("R1")
            r2 = reads.get("r2") or reads.get("R2")
            if r1 and r2:
                return {"r1": r1, "r2": r2}
        if isinstance(reads, list) and len(reads) >= 2:
            return {"r1": reads[0], "r2": reads[1]}
        raise ValueError("Illumina paired-end read_set requires reads.r1/r2 (or R1/R2), or a two-file list.")
    if rtype == "illumina_se":
        if isinstance(reads, str):
            return [reads]
        if isinstance(reads, list) and len(reads) >= 1:
            return reads
        raise ValueError("Illumina single-end read_set requires a read path string or list with one file.")
    if isinstance(reads, str):
        return [reads]
    if isinstance(reads, list):
        return reads
    raise ValueError(f"Unsupported reads format for read_set type={rtype}.")


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
            if reads is None and r.get("source"):
                reads = get_by_dotted_path(config, str(r["source"]))
                if reads is None:
                    raise ValueError(f"read_support.read_sets[{i}] source not found: {r['source']}")
            if reads is None:
                continue
            reads = _normalize_reads_for_type(reads, rtype)
            out.append({"name": name, "type": rtype, "mapper": mapper, "preset": preset, "reads": reads})
        return out

    # new centralized top-level read declarations
    out = []
    for grp in ("long", "short"):
        entries = safe_get(config, ["reads", grp], None) or []
        if not isinstance(entries, list):
            continue
        for i, r in enumerate(entries):
            if not isinstance(r, dict):
                continue
            name = r.get("name") or f"{grp}{i+1}"
            rtype = (r.get("type") or ("pacbio_hifi" if grp == "long" else "illumina_pe")).lower()
            mapper = r.get("mapper", "minimap2")
            preset = r.get("preset") or DEFAULT_PRESET.get(rtype, "map-hifi")
            reads = r.get("reads")
            if reads is None and grp == "short":
                reads = {"r1": r.get("r1") or r.get("R1"), "r2": r.get("r2") or r.get("R2")}
            if reads is None:
                continue
            reads = _normalize_reads_for_type(reads, rtype)
            out.append({"name": name, "type": rtype, "mapper": mapper, "preset": preset, "reads": reads})
    if out:
        return out

    # backward compatibility
    out = []
    if bool(safe_get(config, ["read_support", "use_hifi"], False)):
        hifi = safe_get(config, ["read_support", "hifi_reads"], []) or []
        if hifi:
            out.append({"name": "hifi", "type": "pacbio_hifi", "mapper": "minimap2", "preset": "map-hifi", "reads": hifi})
    if bool(safe_get(config, ["read_support", "use_illumina"], False)):
        il = safe_get(config, ["read_support", "illumina_reads"], {}) or {}
        r1 = il.get("r1") or il.get("R1")
        r2 = il.get("r2") or il.get("R2")
        if r1 and r2:
            out.append({"name": "illumina", "type": "illumina_pe", "mapper": "minimap2", "preset": "sr", "reads": {"r1": r1, "r2": r2}})
    return out


def generate_readset_consensus_recommendations(config: dict, read_support_dir: Path):
    in_tsv = read_support_dir / "problematic_stop_read_support.tsv"
    out_tsv = read_support_dir / "readset_consensus_recommendations.tsv"
    out_md = read_support_dir / "readset_consensus_recommendations.md"
    cols = [
        "gene","stop_aa_position","genomic_codon_start","genomic_codon_end","reference_codon",
        "read_sets_evaluated","read_sets_stop_supported","read_sets_possible_error","read_sets_low_coverage","read_sets_ambiguous",
        "major_codon_by_read_set","major_aa_by_read_set","recommendation_by_read_set","consensus_recommendation","priority","evidence_summary","comment"
    ]
    if not in_tsv.exists():
        pd.DataFrame(columns=cols).to_csv(out_tsv, sep="\t", index=False)
        out_md.write_text("# Read-set consensus recommendations\n\nNo read support TSV found.\n", encoding="utf-8")
        return
    df = pd.read_csv(in_tsv, sep="\t")
    if df.empty:
        pd.DataFrame(columns=cols).to_csv(out_tsv, sep="\t", index=False)
        out_md.write_text("# Read-set consensus recommendations\n\nNo rows available.\n", encoding="utf-8")
        return
    long_types = {"pacbio_hifi","hifi","pacbio_clr","clr","ont","nanopore"}
    type_map = {r["name"]: r["type"] for r in resolve_read_sets(config)}
    rows = []
    grp_cols = ["gene","stop_aa_position","genomic_codon_start","genomic_codon_end","reference_codon"]
    for key, g in df.groupby(grp_cols, dropna=False):
        n = len(g)
        stop_sup = int((g["recommendation"] == "STOP_SUPPORTED_BY_READS").sum())
        poss_err = int((g["recommendation"] == "POSSIBLE_POLISHING_OR_ASSEMBLY_ERROR").sum())
        low_cov = int((g["recommendation"] == "LOW_COVERAGE_REVIEW").sum())
        amb = int((g["recommendation"] == "AMBIGUOUS_READ_SUPPORT").sum())
        has_stop = stop_sup > 0
        has_err = poss_err > 0
        long_err = any((type_map.get(r.read_set, "").lower() in long_types) and (r.recommendation == "POSSIBLE_POLISHING_OR_ASSEMBLY_ERROR") for r in g.itertuples())
        short_err = any((type_map.get(r.read_set, "").lower() not in long_types) and (r.recommendation == "POSSIBLE_POLISHING_OR_ASSEMBLY_ERROR") for r in g.itertuples())
        long_stop = any((type_map.get(r.read_set, "").lower() in long_types) and (r.recommendation == "STOP_SUPPORTED_BY_READS") for r in g.itertuples())
        short_stop = any((type_map.get(r.read_set, "").lower() not in long_types) and (r.recommendation == "STOP_SUPPORTED_BY_READS") for r in g.itertuples())
        if n > 0 and stop_sup == n:
            consensus, priority = "STOP_CONFIRMED_BY_ALL_READSETS", "HIGH"
        elif n > 0 and poss_err == n:
            consensus, priority = "CORRECTION_SUPPORTED_BY_ALL_READSETS", "HIGH"
        elif has_stop and has_err:
            consensus, priority = "CONFLICTING_READSET_EVIDENCE", "HIGH"
        elif has_stop and (stop_sup < n) and (not has_err):
            consensus, priority = "STOP_SUPPORTED_BY_SOME_READSETS", "MEDIUM"
        elif has_err and (poss_err < n) and (not has_stop):
            consensus, priority = "CORRECTION_SUPPORTED_BY_SOME_READSETS", "MEDIUM"
        elif long_err and short_stop:
            consensus, priority = "LONG_READ_ONLY_CORRECTION_SIGNAL", "MEDIUM"
        elif short_err and long_stop:
            consensus, priority = "SHORT_READ_ONLY_CORRECTION_SIGNAL", "MEDIUM"
        elif low_cov > 0 and not (has_stop or has_err):
            consensus, priority = "INSUFFICIENT_READ_SUPPORT", "MEDIUM"
        elif amb > 0:
            consensus, priority = "AMBIGUOUS_EVIDENCE_MANUAL_REVIEW", "MEDIUM"
        else:
            consensus, priority = "MANUAL_REVIEW", "MEDIUM"
        by_set = lambda c: ";".join([f"{r.read_set}:{getattr(r,c)}" for r in g.itertuples()])
        rows.append({
            "gene": key[0], "stop_aa_position": key[1], "genomic_codon_start": key[2], "genomic_codon_end": key[3], "reference_codon": key[4],
            "read_sets_evaluated": n, "read_sets_stop_supported": stop_sup, "read_sets_possible_error": poss_err, "read_sets_low_coverage": low_cov, "read_sets_ambiguous": amb,
            "major_codon_by_read_set": by_set("major_read_codon"), "major_aa_by_read_set": by_set("major_read_amino_acid"),
            "recommendation_by_read_set": by_set("recommendation"), "consensus_recommendation": consensus, "priority": priority,
            "evidence_summary": f"stop_supported={stop_sup}; possible_error={poss_err}; low_coverage={low_cov}; ambiguous={amb}",
            "comment": "diagnostic-only consensus across read sets"
        })
    out_df = pd.DataFrame(rows, columns=cols)
    out_df.to_csv(out_tsv, sep="\t", index=False)
    with open(out_md, "w", encoding="utf-8") as md:
        md.write("# Read-set consensus recommendations\n\n")
        md.write(f"- Total stop sites avaliados: {len(out_df)}\n\n")
        md.write("## Tabela resumida\n\n")
        md.write("|gene|stop_pos|consensus|priority|\n|---|---:|---|---|\n")
        for r in out_df.itertuples():
            md.write(f"|{r.gene}|{r.stop_aa_position}|{r.consensus_recommendation}|{r.priority}|\n")
        md.write("\n## Conflitos entre read sets\n\n")
        conflicts = out_df[out_df["consensus_recommendation"] == "CONFLICTING_READSET_EVIDENCE"]
        md.write(f"- Sites com conflito: {len(conflicts)}\n")
        md.write("\n## Casos insuficientes/ambíguos\n\n")
        ins = out_df[out_df["consensus_recommendation"] == "INSUFFICIENT_READ_SUPPORT"]
        ambg = out_df[out_df["consensus_recommendation"] == "AMBIGUOUS_EVIDENCE_MANUAL_REVIEW"]
        md.write(f"- INSUFFICIENT_READ_SUPPORT: {len(ins)}\n")
        md.write(f"- AMBIGUOUS_EVIDENCE_MANUAL_REVIEW: {len(ambg)}\n")
        md.write("\n## Correções apoiadas por todos\n\n")
        all_corr = out_df[out_df["consensus_recommendation"] == "CORRECTION_SUPPORTED_BY_ALL_READSETS"]
        md.write(f"- Sites: {len(all_corr)}\n")
        md.write("\n## Stops apoiados por alguns read sets\n\n")
        some_stop = out_df[out_df["consensus_recommendation"] == "STOP_SUPPORTED_BY_SOME_READSETS"]
        md.write(f"- Sites: {len(some_stop)}\n")
        md.write("\n## Correções apoiadas por alguns read sets\n\n")
        some_corr = out_df[out_df["consensus_recommendation"] == "CORRECTION_SUPPORTED_BY_SOME_READSETS"]
        md.write(f"- Sites: {len(some_corr)}\n")
        md.write("\n## Stops confirmados por todos\n\n")
        all_stop = out_df[out_df["consensus_recommendation"] == "STOP_CONFIRMED_BY_ALL_READSETS"]
        md.write(f"- Sites: {len(all_stop)}\n\n")
        md.write("> Nenhuma alteração automática foi aplicada ao GenBank.\n")


def _extract_ref_fasta(refined_gb: Path, out_fa: Path):
    rec, _ = read_record(refined_gb)
    SeqIO.write(rec, str(out_fa), "fasta")


def _build_bam_for_read_set(read_set: dict, refined_fa: Path, outdir: Path, threads: int):
    name, preset = read_set["name"], read_set["preset"]
    bam = outdir / f"{name}_to_refined.bam"
    bai = Path(f"{bam}.bai")
    reuse_existing = bool(read_set.get("reuse_existing_bam", True))
    force_remap = bool(read_set.get("force_remap", False))
    if reuse_existing and (not force_remap) and bam.exists() and bai.exists():
        print(f"Reusing existing BAM for read_set {name}: {bam}")
        return bam
    reads = read_set["reads"]
    if isinstance(reads, dict):
        read_args = f"{reads.get('r1','')} {reads.get('r2','')}".strip()
    else:
        read_args = " ".join(reads)
    print(f"Mapping read_set {name} with minimap2 preset {preset} ...")
    cmd = f"minimap2 -t {threads} -ax {preset} {refined_fa} {read_args} | samtools sort -o {bam}"
    run_cmd(["bash", "-lc", cmd], check=False)
    run_cmd(["samtools", "index", str(bam)], check=False)
    return bam if bam.exists() else None


def _aa_for_codon(codon: str, code: int):
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "X"
    return str(Seq(codon).translate(table=code, to_stop=False))


def reverse_complement_dna(seq: str) -> str:
    t = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(t)[::-1].upper()


def run_read_support(config: dict, refined_gb: Path, refinement_dir: Path, outdir: Path):
    import pysam  # lazy import

    outdir = ensure_dir(outdir)
    stop_tsv = refinement_dir / "problematic_cds_stop_context.tsv"
    support_tsv = outdir / "problematic_stop_read_support.tsv"
    var_tsv = outdir / "problematic_stop_variants.tsv"
    summary_md = outdir / "read_support_summary.md"

    sup_cols = ["read_set","gene","seqid","cds_start","cds_end","strand","stop_aa_position","genomic_codon_start","genomic_codon_end","reference_codon","translation_table","read_depth_min","reads_covering_full_codon","major_read_codon","major_read_codon_count","major_read_codon_frequency","major_read_amino_acid","would_remove_stop","alternative_codons","alternative_codon_counts","recommendation","comment"]
    var_cols = ["read_set","gene","seqid","cds_strand","genomic_position","reference_base","depth","A","C","G","T","N","major_base","major_base_frequency","minor_bases","comment"]

    if not stop_tsv.exists():
        pd.DataFrame(columns=sup_cols).to_csv(support_tsv, sep="\t", index=False)
        pd.DataFrame(columns=var_cols).to_csv(var_tsv, sep="\t", index=False)
        summary_md.write_text("# Read support summary\n\nNo problematic stop context file found.\n", encoding="utf-8")
        generate_readset_consensus_recommendations(config, outdir)
        return outdir

    read_sets = resolve_read_sets(config)
    reuse_existing = bool(safe_get(config, ["read_support", "reuse_existing_bam"], True))
    force_remap = bool(safe_get(config, ["read_support", "force_remap"], False))
    read_sets = [{**rs, "reuse_existing_bam": reuse_existing, "force_remap": force_remap} for rs in read_sets]
    refined_fa = outdir / "refined.fa"
    _extract_ref_fasta(refined_gb, refined_fa)
    threads = int(safe_get(config, ["read_support", "threads"], 8))
    bams = {rs["name"]: _build_bam_for_read_set(rs, refined_fa, outdir, threads) for rs in read_sets}

    code = get_genetic_code(config, default=5)
    flank = int(safe_get(config, ["read_support", "flank_bp"], 20))
    stops = pd.read_csv(stop_tsv, sep="\t")
    support_rows, var_rows = [], []

    ref_record, _ = read_record(refined_gb)
    ref_seq = str(ref_record.seq).upper()

    for rs in read_sets:
        bam = bams.get(rs["name"])
        aln = pysam.AlignmentFile(str(bam), "rb") if bam and bam.exists() else None
        for _, r in stops.iterrows():
            g1, g2 = int(r["genomic_codon_start"]), int(r["genomic_codon_end"])
            strand = str(r.get("strand", "+"))
            ref_codon = str(r.get("codon", "NNN"))
            codon_counter, depth_per_pos = Counter(), []
            ref_positions = [g1 - 1, g1, g1 + 1]  # 0-based codon positions

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
                    ref_base = ref_seq[p - 1] if 0 <= (p - 1) < len(ref_seq) else "N"
                    var_rows.append({"read_set": rs["name"], "gene": r["gene"], "seqid": r["seqid"], "cds_strand": strand, "genomic_position": p, "reference_base": ref_base, "depth": depth, "A": counts.get("A",0), "C": counts.get("C",0), "G": counts.get("G",0), "T": counts.get("T",0), "N": counts.get("N",0), "major_base": major, "major_base_frequency": round(freq,3), "minor_bases": minors if minors else ".", "comment": "pileup base counts in genomic orientation"})

                for read in aln.fetch(str(r["seqid"]), max(0, g1 - 1 - flank), g2 + flank):
                    if not read.query_sequence:
                        continue
                    base_by_ref_pos = {}
                    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                        if ref_pos not in ref_positions:
                            continue
                        if query_pos is None:
                            continue
                        if query_pos < 0 or query_pos >= len(read.query_sequence):
                            continue
                        b = read.query_sequence[query_pos].upper()
                        if b not in "ACGT":
                            b = "N"
                        base_by_ref_pos[ref_pos] = b
                    if all(pos in base_by_ref_pos for pos in ref_positions):
                        codon_genomic = "".join(base_by_ref_pos[pos] for pos in sorted(ref_positions))
                        codon_cds = codon_genomic if strand != "-" else reverse_complement_dna(codon_genomic)
                        codon_counter[codon_cds] += 1

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
    generate_readset_consensus_recommendations(config, outdir)
    return outdir
