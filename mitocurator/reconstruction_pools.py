from __future__ import annotations
from pathlib import Path
import gzip
import random
import pandas as pd

from .utils import ensure_dir, safe_get
from .read_support import resolve_read_sets

LONG_TYPES = {"pacbio_hifi", "hifi", "pacbio_clr", "clr", "ont", "nanopore", "illumina_se", "se"}
PE_TYPES = {"illumina_pe", "illumina", "pe"}


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
            seq = fh.readline().rstrip(); fh.readline(); qual = fh.readline().rstrip()
            reads.append((h[1:].split()[0], seq, qual))
    return reads


def _resolve_readset_types(config):
    m = {}
    for rs in resolve_read_sets(config):
        m[str(rs["name"])] = str(rs.get("type", "")).lower()
    # from reads.long/short
    for grp in (safe_get(config, ["reads", "long"], []) or []):
        if isinstance(grp, dict):
            m[str(grp.get("name", ""))] = str(grp.get("type", "")).lower()
    for grp in (safe_get(config, ["reads", "short"], []) or []):
        if isinstance(grp, dict):
            m[str(grp.get("name", ""))] = str(grp.get("type", "")).lower()
    return m


def _infer_output_format_from_readset(read_set: str, type_map: dict):
    rt = type_map.get(read_set, "").lower()
    if rt in PE_TYPES:
        return "paired_fastq", ""
    if rt in LONG_TYPES:
        return "single_fastq", ""
    n = read_set.lower()
    if any(k in n for k in ["illumina", "_pe", "-pe", " pe"]):
        return "paired_fastq", ""
    if any(k in n for k in ["hifi", "pacbio", "clr", "ont", "nanopore"]):
        return "single_fastq", ""
    return "single_fastq", "warning: ambiguous read_set type, defaulting to single_fastq"


def run_reconstruction_pools(config, root: Path, read_support_dir: Path, targeted_extraction_dir: Path, outdir: Path):
    import pysam  # lazy import
    outdir = ensure_dir(outdir)
    tsv_in = targeted_extraction_dir / "targeted_read_extraction.tsv"
    tsv_out = outdir / "reconstruction_pools.tsv"
    md_out = outdir / "reconstruction_pools.md"
    cols = ["target_id", "target_type", "gene", "read_set", "pool_type", "input_sources", "output_fastq", "output_fastq_r1", "output_fastq_r2", "output_format", "reads_or_pairs_written", "singletons_written", "downsampled", "status", "comment"]
    if not tsv_in.exists():
        pd.DataFrame(columns=cols).to_csv(tsv_out, sep="\t", index=False)
        md_out.write_text("# Reconstruction pools\n\nNo targeted extraction TSV found.\n", encoding="utf-8")
        return outdir
    df = pd.read_csv(tsv_in, sep="\t").fillna(".")
    type_map = _resolve_readset_types(config)

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
        tdir = ensure_dir(outdir / target_id / "reads")
        src_fmt = str(r.get("output_format", ".")).strip()
        inferred_fmt, warn = _infer_output_format_from_readset(read_set, type_map)
        fmt = src_fmt if src_fmt in {"single_fastq", "paired_fastq", "interleaved_fastq"} else inferred_fmt

        # load target_only source reads from 08
        ofq, ofq1, ofq2 = Path(str(r.get("output_fastq", "."))), Path(str(r.get("output_fastq_r1", "."))), Path(str(r.get("output_fastq_r2", ".")))
        t_single, t_r1, t_r2 = [], [], []
        if fmt == "paired_fastq" and str(ofq1) != "." and str(ofq2) != ".":
            t_r1, t_r2 = _read_fastq(ofq1), _read_fastq(ofq2)
        elif str(ofq) != ".":
            t_single = _read_fastq(ofq)

        if make_target_only:
            if fmt == "paired_fastq":
                out_r1 = tdir / f"{target_id}.{read_set}.target_only_R1.fastq.gz"
                out_r2 = tdir / f"{target_id}.{read_set}.target_only_R2.fastq.gz"
                if not (reuse and out_r1.exists() and out_r2.exists()):
                    _write_fastq_pe(out_r1, out_r2, t_r1, t_r2)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "target_only", "input_sources": "08_targeted_extraction", "output_fastq": ".", "output_fastq_r1": str(out_r1), "output_fastq_r2": str(out_r2), "output_format": "paired_fastq", "reads_or_pairs_written": min(len(t_r1), len(t_r2)), "singletons_written": abs(len(t_r1)-len(t_r2)), "downsampled": "no", "status": "ok", "comment": "from target extraction R1/R2" + ("; " + warn if warn else "")})
            elif fmt == "interleaved_fastq":
                out_fq = tdir / f"{target_id}.{read_set}.target_only.interleaved.fastq.gz"
                if not (reuse and out_fq.exists()):
                    _write_fastq(out_fq, t_single)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "target_only", "input_sources": "08_targeted_extraction", "output_fastq": str(out_fq), "output_fastq_r1": ".", "output_fastq_r2": ".", "output_format": "interleaved_fastq", "reads_or_pairs_written": len(t_single), "singletons_written": len(t_single), "downsampled": "no", "status": "ok", "comment": "from interleaved target extraction" + ("; " + warn if warn else "")})
            else:
                out_fq = tdir / f"{target_id}.{read_set}.target_only.fastq.gz"
                if not (reuse and out_fq.exists()):
                    _write_fastq(out_fq, t_single)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "target_only", "input_sources": "08_targeted_extraction", "output_fastq": str(out_fq), "output_fastq_r1": ".", "output_fastq_r2": ".", "output_format": "single_fastq", "reads_or_pairs_written": len(t_single), "singletons_written": len(t_single), "downsampled": "no", "status": "ok", "comment": "from target extraction single" + ("; " + warn if warn else "")})

        bam = read_support_dir / f"{read_set}_to_refined.bam"
        bai = Path(str(bam) + ".bai")
        m_single, m_r1, m_r2, m_singletons = [], [], [], []
        down = "no"
        if (make_mito or make_combined) and bam.exists() and bai.exists():
            aln = pysam.AlignmentFile(str(bam), "rb")
            if fmt == "paired_fastq":
                pairs = {}
                for read in aln.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary or int(read.mapping_quality) < min_mapq:
                        continue
                    if not read.query_sequence:
                        continue
                    d = pairs.setdefault(read.query_name, {})
                    q = read.qual if read.qual else ("I" * len(read.query_sequence))
                    if read.is_read1:
                        d["r1"] = (read.query_name, read.query_sequence, q)
                    elif read.is_read2:
                        d["r2"] = (read.query_name, read.query_sequence, q)
                names = list(pairs.keys())
                if len(names) > max_short:
                    names = rng.sample(names, max_short); down = "yes"
                for n in names:
                    d = pairs[n]
                    if "r1" in d and "r2" in d:
                        m_r1.append(d["r1"]); m_r2.append(d["r2"])
                    elif "r1" in d:
                        m_singletons.append(d["r1"])
                    elif "r2" in d:
                        m_singletons.append(d["r2"])
            else:
                seen = set()
                for read in aln.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary or int(read.mapping_quality) < min_mapq:
                        continue
                    if not read.query_sequence or read.query_name in seen:
                        continue
                    seen.add(read.query_name)
                    q = read.qual if read.qual else ("I" * len(read.query_sequence))
                    m_single.append((read.query_name, read.query_sequence, q))
                if len(m_single) > max_long:
                    m_single = rng.sample(m_single, max_long); down = "yes"
            aln.close()

        if make_mito:
            if fmt == "paired_fastq":
                out_r1 = tdir / f"{target_id}.{read_set}.mitogenome_mapped_R1.fastq.gz"
                out_r2 = tdir / f"{target_id}.{read_set}.mitogenome_mapped_R2.fastq.gz"
                if not (reuse and out_r1.exists() and out_r2.exists()):
                    _write_fastq_pe(out_r1, out_r2, m_r1, m_r2)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "mitogenome_mapped", "input_sources": "06_read_support BAM", "output_fastq": ".", "output_fastq_r1": str(out_r1), "output_fastq_r2": str(out_r2), "output_format": "paired_fastq", "reads_or_pairs_written": min(len(m_r1), len(m_r2)), "singletons_written": len(m_singletons), "downsampled": down, "status": "ok" if bam.exists() else "missing_bam", "comment": "mapped reads from whole mitogenome"})
            else:
                out_fq = tdir / f"{target_id}.{read_set}.mitogenome_mapped.fastq.gz"
                if not (reuse and out_fq.exists()):
                    _write_fastq(out_fq, m_single)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "mitogenome_mapped", "input_sources": "06_read_support BAM", "output_fastq": str(out_fq), "output_fastq_r1": ".", "output_fastq_r2": ".", "output_format": "single_fastq", "reads_or_pairs_written": len(m_single), "singletons_written": len(m_single), "downsampled": down, "status": "ok" if bam.exists() else "missing_bam", "comment": "mapped reads from whole mitogenome"})

        if make_combined:
            if fmt == "paired_fastq":
                by = {}
                for x in (t_r1 + m_r1): by.setdefault(x[0], {})["r1"] = x
                for x in (t_r2 + m_r2): by.setdefault(x[0], {})["r2"] = x
                c_r1, c_r2, c_single = [], [], []
                for n, d in by.items():
                    if "r1" in d and "r2" in d:
                        c_r1.append(d["r1"]); c_r2.append(d["r2"])
                    else:
                        c_single.extend(list(d.values()))
                if len(c_r1) > max_short:
                    idx = set(rng.sample(range(len(c_r1)), max_short))
                    c_r1 = [x for i, x in enumerate(c_r1) if i in idx]
                    c_r2 = [x for i, x in enumerate(c_r2) if i in idx]
                    downc = "yes"
                else:
                    downc = "no"
                out_r1 = tdir / f"{target_id}.{read_set}.combined_R1.fastq.gz"
                out_r2 = tdir / f"{target_id}.{read_set}.combined_R2.fastq.gz"
                if not (reuse and out_r1.exists() and out_r2.exists()):
                    _write_fastq_pe(out_r1, out_r2, c_r1, c_r2)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "combined", "input_sources": "target_only+mitogenome_mapped", "output_fastq": ".", "output_fastq_r1": str(out_r1), "output_fastq_r2": str(out_r2), "output_format": "paired_fastq", "reads_or_pairs_written": min(len(c_r1), len(c_r2)), "singletons_written": len(c_single), "downsampled": downc, "status": "ok", "comment": "deduplicated union"})
            else:
                by = {x[0]: x for x in (t_single + m_single)}
                c = list(by.values())
                if len(c) > max_long:
                    c = rng.sample(c, max_long); downc = "yes"
                else:
                    downc = "no"
                out_fq = tdir / f"{target_id}.{read_set}.combined.fastq.gz"
                if not (reuse and out_fq.exists()):
                    _write_fastq(out_fq, c)
                rows.append({"target_id": target_id, "target_type": r["target_type"], "gene": r["gene"], "read_set": read_set, "pool_type": "combined", "input_sources": "target_only+mitogenome_mapped", "output_fastq": str(out_fq), "output_fastq_r1": ".", "output_fastq_r2": ".", "output_format": "single_fastq", "reads_or_pairs_written": len(c), "singletons_written": len(c), "downsampled": downc, "status": "ok", "comment": "deduplicated union"})

    out_df = pd.DataFrame(rows, columns=cols)
    out_df.to_csv(tsv_out, sep="\t", index=False)
    with open(md_out, "w", encoding="utf-8") as md:
        md.write("# Reconstruction pools\n\n")
        md.write(f"- total de targets: {out_df['target_id'].nunique() if not out_df.empty else 0}\n")
        md.write(f"- total de pools criados: {len(out_df)}\n\n")
        if not out_df.empty:
            md.write(f"- single_fastq pools: {int((out_df['output_format']=='single_fastq').sum())}\n")
            md.write(f"- paired_fastq pools: {int((out_df['output_format']=='paired_fastq').sum())}\n")
            md.write(f"- interleaved_fastq pools: {int((out_df['output_format']=='interleaved_fastq').sum())}\n\n")
        md.write("## Resumo por pool_type\n\n")
        for k, g in (out_df.groupby("pool_type") if not out_df.empty else []):
            md.write(f"- {k}: {len(g)}\n")
        md.write("\n## Resumo por read_set\n\n")
        for k, g in (out_df.groupby("read_set") if not out_df.empty else []):
            md.write(f"- {k}: {len(g)} pools\n")
        md.write("\n## Alvos sem reads\n\n")
        nr = out_df[out_df["reads_or_pairs_written"] == 0] if not out_df.empty else pd.DataFrame()
        if nr.empty:
            md.write("- Nenhum\n")
        else:
            for rr in nr.itertuples():
                md.write(f"- {rr.target_id} [{rr.read_set}/{rr.pool_type}]\n")
        md.write("\n> Nenhum montador foi executado e nenhuma correção automática foi aplicada ao GenBank.\n")
    return outdir
