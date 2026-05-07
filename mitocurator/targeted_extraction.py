from __future__ import annotations
from pathlib import Path
import gzip
import random
import pandas as pd

from .utils import ensure_dir, safe_get
from .read_support import resolve_read_sets


def _collect_targets(config: dict, refinement_dir: Path, read_support_dir: Path):
    flank = int(safe_get(config, ["targeted_extraction", "flank_bp"], 3000))
    include_missing = bool(safe_get(config, ["targeted_extraction", "include_missing_gene_candidates"], True))
    include_problematic = bool(safe_get(config, ["targeted_extraction", "include_problematic_cds"], True))
    include_consensus = bool(safe_get(config, ["targeted_extraction", "include_consensus_correction_candidates"], True))
    targets = []

    if include_missing:
        tsv = refinement_dir / "reference_similarity_candidates.tsv"
        if tsv.exists():
            df = pd.read_csv(tsv, sep="\t")
            sub = df[df["decision_hint"].isin(["PARTIAL_REF_MATCH", "STRONG_REF_MATCH"])]
            for _, r in sub.iterrows():
                targets.append({
                    "seqid": str(r["seqid"]),
                    "start1": int(r["start"]),
                    "end1": int(r["end"]),
                    "target_id": f"missing_{r['gene']}_{r['candidate_id']}",
                    "target_type": "missing_gene_candidate",
                    "gene": str(r["gene"]),
                    "comment": f"from reference_similarity_candidates ({r['decision_hint']})",
                    "flank_bp": flank,
                })

    if include_problematic:
        tsv = refinement_dir / "problematic_cds_stop_context.tsv"
        if tsv.exists():
            df = pd.read_csv(tsv, sep="\t")
            for gene, g in df.groupby("gene"):
                first = g.iloc[0]
                targets.append({
                    "seqid": str(first["seqid"]),
                    "start1": int(first["cds_start"]),
                    "end1": int(first["cds_end"]),
                    "target_id": f"problematic_{gene}",
                    "target_type": "problematic_cds",
                    "gene": str(gene),
                    "comment": "CDS with internal stop(s)",
                    "flank_bp": flank,
                })

    if include_consensus:
        tsv = read_support_dir / "readset_consensus_recommendations.tsv"
        if tsv.exists():
            df = pd.read_csv(tsv, sep="\t")
            keep = {"CORRECTION_SUPPORTED_BY_ALL_READSETS", "CORRECTION_SUPPORTED_BY_SOME_READSETS", "CONFLICTING_READSET_EVIDENCE"}
            sub = df[df["consensus_recommendation"].isin(keep)]
            # infer seqid from problematic context when possible
            ctx_tsv = refinement_dir / "problematic_cds_stop_context.tsv"
            seqid_map = {}
            if ctx_tsv.exists():
                ctx = pd.read_csv(ctx_tsv, sep="\t")
                for _, rr in ctx.iterrows():
                    seqid_map[(str(rr["gene"]), int(rr["stop_aa_position"]))] = str(rr["seqid"])
            for _, r in sub.iterrows():
                s1 = int(r["genomic_codon_start"])
                e1 = int(r["genomic_codon_end"])
                gene = str(r["gene"])
                stop = int(r["stop_aa_position"])
                targets.append({
                    "seqid": seqid_map.get((gene, stop), "."),
                    "start1": s1,
                    "end1": e1,
                    "target_id": f"correction_{gene}_stop{stop}",
                    "target_type": "consensus_correction_candidate",
                    "gene": gene,
                    "comment": str(r["consensus_recommendation"]),
                    "flank_bp": flank,
                })
    return targets


def _write_fastq_gz(path: Path, reads):
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for name, seq, qual in reads:
            out.write(f"@{name}\n{seq}\n+\n{qual}\n")


def run_targeted_extraction(config, root: Path, refinement_dir: Path, read_support_dir: Path, outdir: Path):
    import pysam  # lazy import
    outdir = ensure_dir(outdir)
    reads_dir = ensure_dir(outdir / "reads")
    targets_bed = outdir / "targets.bed"
    out_tsv = outdir / "targeted_read_extraction.tsv"
    out_md = outdir / "targeted_read_extraction.md"
    cols = ["target_id","target_type","gene","seqid","start","end","flank_bp","read_set","bam","reads_extracted","output_fastq","status","comment"]

    reuse_outputs = bool(safe_get(config, ["targeted_extraction", "reuse_existing_outputs"], True))
    max_reads = int(safe_get(config, ["targeted_extraction", "max_reads_per_target"], 5000))
    min_mapq = int(safe_get(config, ["targeted_extraction", "min_mapping_quality"], 20))
    seed = int(safe_get(config, ["targeted_extraction", "random_seed"], 42))
    rng = random.Random(seed)

    targets = _collect_targets(config, refinement_dir, read_support_dir)
    if not targets:
        pd.DataFrame(columns=cols).to_csv(out_tsv, sep="\t", index=False)
        targets_bed.write_text("seqid\tstart0\tend\ttarget_id\ttarget_type\tgene\tcomment\n", encoding="utf-8")
        out_md.write_text("# Targeted read extraction\n\nNo targets selected.\n", encoding="utf-8")
        return outdir

    # write bed
    with open(targets_bed, "w", encoding="utf-8") as bed:
        bed.write("seqid\tstart0\tend\ttarget_id\ttarget_type\tgene\tcomment\n")
        for t in targets:
            s0 = max(0, int(t["start1"]) - 1 - int(t["flank_bp"]))
            e0 = max(int(t["end1"]), int(t["start1"])) + int(t["flank_bp"])
            bed.write(f"{t['seqid']}\t{s0}\t{e0}\t{t['target_id']}\t{t['target_type']}\t{t['gene']}\t{t['comment']}\n")

    rows = []
    read_sets = resolve_read_sets(config)
    for rs in read_sets:
        bam = read_support_dir / f"{rs['name']}_to_refined.bam"
        bai = Path(f"{bam}.bai")
        if not (bam.exists() and bai.exists()):
            for t in targets:
                rows.append({"target_id": t["target_id"], "target_type": t["target_type"], "gene": t["gene"], "seqid": t["seqid"], "start": t["start1"], "end": t["end1"], "flank_bp": t["flank_bp"], "read_set": rs["name"], "bam": str(bam), "reads_extracted": 0, "output_fastq": ".", "status": "missing_bam", "comment": "BAM/BAI not found in 06_read_support"})
            continue
        aln = pysam.AlignmentFile(str(bam), "rb")
        for t in targets:
            s0 = max(0, int(t["start1"]) - 1 - int(t["flank_bp"]))
            e0 = max(int(t["end1"]), int(t["start1"])) + int(t["flank_bp"])
            is_pe = rs["type"] in {"illumina_pe", "illumina", "pe"}
            suffix = f"{t['target_id']}.{rs['name']}.interleaved.fastq.gz" if is_pe else f"{t['target_id']}.{rs['name']}.fastq.gz"
            out_fastq = reads_dir / suffix
            if reuse_outputs and out_fastq.exists():
                rows.append({"target_id": t["target_id"], "target_type": t["target_type"], "gene": t["gene"], "seqid": t["seqid"], "start": t["start1"], "end": t["end1"], "flank_bp": t["flank_bp"], "read_set": rs["name"], "bam": str(bam), "reads_extracted": -1, "output_fastq": str(out_fastq), "status": "reused", "comment": "existing output reused"})
                continue
            seen = set()
            recs = []
            for read in aln.fetch(t["seqid"], s0, e0):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if int(read.mapping_quality) < min_mapq:
                    continue
                if read.query_name in seen:
                    continue
                seen.add(read.query_name)
                seq = read.query_sequence or ""
                if not seq:
                    continue
                qual = read.qual if read.qual else ("I" * len(seq))
                recs.append((read.query_name, seq, qual))
            if len(recs) > max_reads:
                recs = rng.sample(recs, max_reads)
                comment = f"downsampled to max_reads_per_target={max_reads}"
            else:
                comment = "."
            _write_fastq_gz(out_fastq, recs)
            rows.append({"target_id": t["target_id"], "target_type": t["target_type"], "gene": t["gene"], "seqid": t["seqid"], "start": t["start1"], "end": t["end1"], "flank_bp": t["flank_bp"], "read_set": rs["name"], "bam": str(bam), "reads_extracted": len(recs), "output_fastq": str(out_fastq), "status": "ok", "comment": comment + ("; interleaved output for PE" if is_pe else "")})
        aln.close()

    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(out_tsv, sep="\t", index=False)
    no_reads = df[df["reads_extracted"] == 0]
    with open(out_md, "w", encoding="utf-8") as md:
        md.write("# Targeted read extraction\n\n")
        md.write(f"- Targets total: {len(set(df['target_id'])) if not df.empty else 0}\n")
        md.write(f"- Rows (target x read_set): {len(df)}\n\n")
        md.write("## Reads extraídas por target/read_set\n\n")
        for r in df.itertuples():
            md.write(f"- {r.target_id} [{r.read_set}]: {r.reads_extracted} ({r.status})\n")
        md.write("\n## Alvos sem reads\n\n")
        if no_reads.empty:
            md.write("- Nenhum\n")
        else:
            for r in no_reads.itertuples():
                md.write(f"- {r.target_id} [{r.read_set}]\n")
        md.write("\n## Próximo passo\n\n")
        md.write("Use os FASTQs alvo para remontagem local/análise manual. Nenhuma correção automática foi aplicada ao GenBank.\n")
    return outdir

