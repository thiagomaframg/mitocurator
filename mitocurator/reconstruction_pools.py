from __future__ import annotations
from pathlib import Path
import gzip
import random
import pandas as pd

from .utils import ensure_dir, safe_get


def _write_fastq(path: Path, reads):
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for name, seq, qual in reads:
            out.write(f"@{name}\n{seq}\n+\n{qual}\n")


def _write_fastq_pe(r1_path: Path, r2_path: Path, r1_reads, r2_reads):
    with gzip.open(r1_path, "wt", encoding="utf-8") as o1, gzip.open(r2_path, "wt", encoding="utf-8") as o2:
        for name, seq, qual in r1_reads:
            o1.write(f"@{name}\n{seq}\n+\n{qual}\n")
        for name, seq, qual in r2_reads:
            o2.write(f"@{name}\n{seq}\n+\n{qual}\n")


def _read_fastq(path: Path):
    reads = []
    if not path.exists():
        return reads
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        while True:
            h = fh.readline().rstrip()
            if not h:
                break
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            reads.append((h[1:].split()[0], seq, qual))
    return reads


def run_reconstruction_pools(config, root: Path, read_support_dir: Path, targeted_extraction_dir: Path, outdir: Path):
    import pysam  # lazy import
    outdir = ensure_dir(outdir)
    tsv_in = targeted_extraction_dir / "targeted_read_extraction.tsv"
    tsv_out = outdir / "reconstruction_pools.tsv"
    md_out = outdir / "reconstruction_pools.md"
    cols = ["target_id","target_type","gene","read_set","pool_type","input_sources","output_fastq","output_fastq_r1","output_fastq_r2","reads_or_pairs_written","singletons_written","downsampled","status","comment"]
    if not tsv_in.exists():
        pd.DataFrame(columns=cols).to_csv(tsv_out, sep="\t", index=False)
        md_out.write_text("# Reconstruction pools\n\nNo targeted extraction TSV found.\n", encoding="utf-8")
        return outdir

    df = pd.read_csv(tsv_in, sep="\t")
    if df.empty:
        pd.DataFrame(columns=cols).to_csv(tsv_out, sep="\t", index=False)
        md_out.write_text("# Reconstruction pools\n\nNo targeted extraction rows found.\n", encoding="utf-8")
        return outdir

    make_target_only = bool(safe_get(config, ["reconstruction_pools", "make_target_only"], True))
    make_mito = bool(safe_get(config, ["reconstruction_pools", "make_mitogenome_mapped"], True))
    make_combined = bool(safe_get(config, ["reconstruction_pools", "make_combined"], True))
    max_long = int(safe_get(config, ["reconstruction_pools", "max_long_reads_per_pool"], 10000))
    max_short = int(safe_get(config, ["reconstruction_pools", "max_short_pairs_per_pool"], 20000))
    min_mapq = int(safe_get(config, ["reconstruction_pools", "min_mapping_quality"], 20))
    reuse = bool(safe_get(config, ["reconstruction_pools", "reuse_existing_outputs"], True))
    rng = random.Random(int(safe_get(config, ["reconstruction_pools", "random_seed"], 42)))

    rows = []
    for _, r in df.iterrows():
        target_id, read_set = str(r["target_id"]), str(r["read_set"])
        target_dir = ensure_dir(outdir / target_id / "reads")
        is_pe = str(r.get("output_format", "")).strip() in {"paired_fastq", "interleaved_fastq"} or read_set.lower().find("illumina") >= 0
        # collect target-only reads from 08 outputs
        t_r1, t_r2, t_single = [], [], []
        ofq = str(r.get("output_fastq", "."))
        ofq1 = str(r.get("output_fastq_r1", "."))
        ofq2 = str(r.get("output_fastq_r2", "."))
        if ofq1 != "." and ofq2 != ".":
            t_r1 = _read_fastq(Path(ofq1))
            t_r2 = _read_fastq(Path(ofq2))
            if ofq != ".":
                t_single = _read_fastq(Path(ofq))
        elif ofq != ".":
            t_single = _read_fastq(Path(ofq))

        if make_target_only:
            if is_pe and t_r1 and t_r2:
                out_r1 = target_dir / f"{target_id}.{read_set}.target_only_R1.fastq.gz"
                out_r2 = target_dir / f"{target_id}.{read_set}.target_only_R2.fastq.gz"
                if not (reuse and out_r1.exists() and out_r2.exists()):
                    _write_fastq_pe(out_r1, out_r2, t_r1, t_r2)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "target_only", "input_sources": "08_targeted_extraction", "output_fastq": ".", "output_fastq_r1": str(out_r1), "output_fastq_r2": str(out_r2), "reads_or_pairs_written": min(len(t_r1), len(t_r2)), "singletons_written": len(t_single), "downsampled": "no", "status": "ok", "comment": "from targeted extraction outputs"})
            else:
                out_fq = target_dir / f"{target_id}.{read_set}.target_only.fastq.gz"
                if not (reuse and out_fq.exists()):
                    _write_fastq(out_fq, t_single)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "target_only", "input_sources": "08_targeted_extraction", "output_fastq": str(out_fq), "output_fastq_r1": ".", "output_fastq_r2": ".", "reads_or_pairs_written": len(t_single), "singletons_written": len(t_single), "downsampled": "no", "status": "ok", "comment": "from targeted extraction outputs"})

        # mitogenome mapped pool from BAM
        m_r1, m_r2, m_single = [], [], []
        if make_mito or make_combined:
            bam = read_support_dir / f"{read_set}_to_refined.bam"
            bai = Path(f"{bam}.bai")
            if bam.exists() and bai.exists():
                aln = pysam.AlignmentFile(str(bam), "rb")
                seen = set()
                pairs = {}
                for read in aln.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if int(read.mapping_quality) < min_mapq:
                        continue
                    seq = read.query_sequence or ""
                    if not seq:
                        continue
                    qual = read.qual if read.qual else ("I" * len(seq))
                    if is_pe:
                        d = pairs.setdefault(read.query_name, {})
                        if read.is_read1:
                            d["r1"] = (read.query_name, seq, qual)
                        elif read.is_read2:
                            d["r2"] = (read.query_name, seq, qual)
                    else:
                        if read.query_name in seen:
                            continue
                        seen.add(read.query_name)
                        m_single.append((read.query_name, seq, qual))
                aln.close()
                if is_pe:
                    names = list(pairs.keys())
                    if len(names) > max_short:
                        names = rng.sample(names, max_short)
                        down = "yes"
                    else:
                        down = "no"
                    for n in names:
                        d = pairs[n]
                        if "r1" in d and "r2" in d:
                            m_r1.append(d["r1"]); m_r2.append(d["r2"])
                        elif "r1" in d:
                            m_single.append(d["r1"])
                        elif "r2" in d:
                            m_single.append(d["r2"])
                else:
                    if len(m_single) > max_long:
                        m_single = rng.sample(m_single, max_long); down = "yes"
                    else:
                        down = "no"
            else:
                down = "no"
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "mitogenome_mapped", "input_sources": "06_read_support BAM", "output_fastq": ".", "output_fastq_r1": ".", "output_fastq_r2": ".", "reads_or_pairs_written": 0, "singletons_written": 0, "downsampled": down, "status": "missing_bam", "comment": "BAM/BAI missing"})
                continue

        if make_mito:
            if is_pe:
                out_r1 = target_dir / f"{target_id}.{read_set}.mitogenome_mapped_R1.fastq.gz"
                out_r2 = target_dir / f"{target_id}.{read_set}.mitogenome_mapped_R2.fastq.gz"
                if not (reuse and out_r1.exists() and out_r2.exists()):
                    _write_fastq_pe(out_r1, out_r2, m_r1, m_r2)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "mitogenome_mapped", "input_sources": "06_read_support BAM", "output_fastq": ".", "output_fastq_r1": str(out_r1), "output_fastq_r2": str(out_r2), "reads_or_pairs_written": min(len(m_r1), len(m_r2)), "singletons_written": len(m_single), "downsampled": down, "status": "ok", "comment": "all mapped reads on refined mitogenome"})
            else:
                out_fq = target_dir / f"{target_id}.{read_set}.mitogenome_mapped.fastq.gz"
                if not (reuse and out_fq.exists()):
                    _write_fastq(out_fq, m_single)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "mitogenome_mapped", "input_sources": "06_read_support BAM", "output_fastq": str(out_fq), "output_fastq_r1": ".", "output_fastq_r2": ".", "reads_or_pairs_written": len(m_single), "singletons_written": len(m_single), "downsampled": down, "status": "ok", "comment": "all mapped reads on refined mitogenome"})

        if make_combined:
            if is_pe:
                by = {}
                for x in t_r1 + m_r1: by.setdefault(x[0], {})["r1"] = x
                for x in t_r2 + m_r2: by.setdefault(x[0], {})["r2"] = x
                c_r1, c_r2, c_single = [], [], []
                for n, d in by.items():
                    if "r1" in d and "r2" in d:
                        c_r1.append(d["r1"]); c_r2.append(d["r2"])
                    else:
                        c_single.extend([d[k] for k in d])
                if len(c_r1) > max_short:
                    idx = set(rng.sample(range(len(c_r1)), max_short))
                    c_r1 = [x for i, x in enumerate(c_r1) if i in idx]
                    c_r2 = [x for i, x in enumerate(c_r2) if i in idx]
                    down_c = "yes"
                else:
                    down_c = "no"
                out_r1 = target_dir / f"{target_id}.{read_set}.combined_R1.fastq.gz"
                out_r2 = target_dir / f"{target_id}.{read_set}.combined_R2.fastq.gz"
                if not (reuse and out_r1.exists() and out_r2.exists()):
                    _write_fastq_pe(out_r1, out_r2, c_r1, c_r2)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "combined", "input_sources": "target_only+mitogenome_mapped", "output_fastq": ".", "output_fastq_r1": str(out_r1), "output_fastq_r2": str(out_r2), "reads_or_pairs_written": min(len(c_r1), len(c_r2)), "singletons_written": len(c_single), "downsampled": down_c, "status": "ok", "comment": "deduplicated union"})
            else:
                by = {x[0]: x for x in (t_single + m_single)}
                c = list(by.values())
                if len(c) > max_long:
                    c = rng.sample(c, max_long); down_c = "yes"
                else:
                    down_c = "no"
                out_fq = target_dir / f"{target_id}.{read_set}.combined.fastq.gz"
                if not (reuse and out_fq.exists()):
                    _write_fastq(out_fq, c)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "combined", "input_sources": "target_only+mitogenome_mapped", "output_fastq": str(out_fq), "output_fastq_r1": ".", "output_fastq_r2": ".", "reads_or_pairs_written": len(c), "singletons_written": len(c), "downsampled": down_c, "status": "ok", "comment": "deduplicated union"})

    out_df = pd.DataFrame(rows, columns=cols)
    out_df.to_csv(tsv_out, sep="\t", index=False)
    with open(md_out, "w", encoding="utf-8") as md:
        md.write("# Reconstruction pools\n\n")
        md.write(f"- total de targets: {out_df['target_id'].nunique() if not out_df.empty else 0}\n")
        md.write(f"- total de pools criados: {len(out_df)}\n\n")
        md.write("## Resumo por pool_type\n\n")
        if not out_df.empty:
            for k, g in out_df.groupby("pool_type"):
                md.write(f"- {k}: {len(g)}\n")
        md.write("\n## Resumo por read_set\n\n")
        if not out_df.empty:
            for k, g in out_df.groupby("read_set"):
                md.write(f"- {k}: {len(g)} pools\n")
        md.write("\n## Alvos sem reads\n\n")
        nr = out_df[out_df["reads_or_pairs_written"] == 0] if not out_df.empty else pd.DataFrame()
        if nr.empty:
            md.write("- Nenhum\n")
        else:
            for r in nr.itertuples():
                md.write(f"- {r.target_id} [{r.read_set}/{r.pool_type}]\n")
        md.write("\n> Nenhum montador foi executado e nenhuma correção automática foi aplicada ao GenBank.\n")
    return outdir

