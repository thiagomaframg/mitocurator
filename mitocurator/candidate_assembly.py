from __future__ import annotations
from pathlib import Path
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .utils import ensure_dir, safe_get, run_cmd, get_genetic_code
from .io import read_record
from .targeted_consensus import _reference_protein

LONG_TYPES = {"hifi","pacbio_hifi","pacbio_clr","clr","ont","nanopore"}
SHORT_TYPES = {"illumina","illumina_pe","pe","illumina_se","se"}


def _type_from_readset(config, read_set):
    for grp in (safe_get(config,["reads","long"],[]) or []) + (safe_get(config,["reads","short"],[]) or []):
        if isinstance(grp, dict) and str(grp.get("name")) == read_set:
            return str(grp.get("type",""))
    for rs in (safe_get(config,["read_support","read_sets"],[]) or []):
        if str(rs.get("name")) == read_set:
            return str(rs.get("type",""))
    n = read_set.lower()
    if "illumina" in n: return "illumina_pe"
    if "hifi" in n: return "pacbio_hifi"
    if "ont" in n: return "ont"
    return ""


def run_candidate_assembly(config, root: Path, consensus_dir: Path, pools_dir: Path, refinement_dir: Path, outdir: Path):
    outdir = ensure_dir(outdir)
    targets_tsv = outdir / "candidate_assembly_targets.tsv"
    sum_tsv = outdir / "candidate_assembly_summary.tsv"
    sum_md = outdir / "candidate_assembly_summary.md"
    ca = safe_get(config,["candidate_assembly"],{}) or {}
    src_mode = ca.get("candidate_source","cross_readset")
    cross = consensus_dir / "cross_readset_missing_gene_candidates.tsv"
    best = consensus_dir / "best_missing_gene_candidates.tsv"
    if src_mode == "cross_readset" and cross.exists():
        cand = pd.read_csv(cross, sep="\t")
        cand = cand.rename(columns={"combined_rank":"rank", "combined_recommendation":"recommendation"})
    elif best.exists():
        cand = pd.read_csv(best, sep="\t")
        cand = cand.rename(columns={"rank":"rank", "recommendation":"recommendation"})
    else:
        pd.DataFrame().to_csv(targets_tsv, sep="\t", index=False)
        pd.DataFrame().to_csv(sum_tsv, sep="\t", index=False)
        sum_md.write_text("# Candidate assembly\n\nNo candidate source found.\n", encoding="utf-8")
        return outdir

    inc = set(ca.get("include_recommendations", ["CROSS_READSET_GENE_CANDIDATE_SUPPORTED","CROSS_READSET_PARTIAL_GENE_CANDIDATE_SUPPORTED","SINGLE_READSET_GENE_CANDIDATE_SUPPORTED"]))
    max_per_gene = int(ca.get("max_candidates_per_gene",1))
    cand = cand[cand["recommendation"].isin(inc)]
    if "rank" in cand.columns:
        cand = cand[cand["rank"] <= max_per_gene]
    cand.to_csv(targets_tsv, sep="\t", index=False)

    pools = pd.read_csv(pools_dir / "reconstruction_pools.tsv", sep="\t") if (pools_dir / "reconstruction_pools.tsv").exists() else pd.DataFrame()
    rec, _ = read_record(refinement_dir / "refined.gb")
    refined_fa = outdir / "refined_mitogenome.fasta"
    SeqIO.write([SeqIO.SeqRecord(rec.seq, id=rec.id, description="")], str(refined_fa), "fasta")
    run_cmd(["makeblastdb", "-in", str(refined_fa), "-dbtype", "nucl"], check=False)
    threads = int(ca.get("threads", 16))
    min_len = int(ca.get("min_contig_len",300))
    code = get_genetic_code(config, default=5)

    rows = []
    for rr in cand.itertuples():
        gene = str(rr.gene); tid = str(rr.target_id)
        base = ensure_dir(outdir / gene / tid)
        ad = ensure_dir(base / "assembly"); bd = ensure_dir(base / "blast"); dd = ensure_dir(base / "diagnosis")
        subset = pools[(pools["target_id"].astype(str)==tid) & (pools["pool_type"].astype(str)=="combined")]
        for pr in subset.itertuples():
            rs = str(pr.read_set); rtype = _type_from_readset(config, rs)
            asm_fa = None; assembler = None; status = "ASSEMBLY_NOT_RUN"
            if rtype in LONG_TYPES and bool(ca.get("run_long_read_assembly", True)):
                assembler = "flye"
                inp = str(pr.output_fastq)
                if inp and inp != ".":
                    mode = "--pacbio-hifi" if "hifi" in rtype else "--pacbio-raw" if "clr" in rtype else "--nano-hq"
                    outd = ad / f"flye_{rs}"
                    run_cmd(["bash","-lc",f"flye {mode} {inp} --out-dir {outd} --threads {threads} --genome-size {safe_get(config,['candidate_assembly','flye','genome_size'],'20k')} {safe_get(config,['candidate_assembly','flye','extra_args'],'')}"], check=False)
                    asm = outd / "assembly.fasta"
                    if asm.exists():
                        asm_fa = ad / f"{tid}.{rs}.flye.fasta"; shutil.copyfile(asm, asm_fa); status = "OK"
                    else:
                        status = "ASSEMBLY_FAILED"
                else:
                    status = "NO_INPUT_READS"
            elif rtype in SHORT_TYPES and bool(ca.get("run_short_read_assembly", True)):
                assembler = "spades"
                outd = ad / f"spades_{rs}"
                if str(pr.output_format) == "paired_fastq" and str(pr.output_fastq_r1)!="." and str(pr.output_fastq_r2)!=".":
                    cmd = f"spades.py -1 {pr.output_fastq_r1} -2 {pr.output_fastq_r2} -o {outd} --threads {threads} {safe_get(config,['candidate_assembly','spades','extra_args'],'')}"
                else:
                    cmd = f"spades.py -s {pr.output_fastq} -o {outd} --threads {threads} {safe_get(config,['candidate_assembly','spades','extra_args'],'')}"
                run_cmd(["bash","-lc",cmd], check=False)
                asm = outd / "contigs.fasta"
                if asm.exists():
                    asm_fa = ad / f"{tid}.{rs}.spades.fasta"; shutil.copyfile(asm, asm_fa); status = "OK"
                else:
                    status = "ASSEMBLY_FAILED"
            else:
                status = "ASSEMBLER_NOT_AVAILABLE"

            best_contig = "."; best_len = 0; bident = bqcov = bsstart = bsend = "."; tbpid=tbqcov=tbscov=tbbs="."
            recm = "NO_LOCAL_ASSEMBLY_SUPPORT"
            if asm_fa and Path(asm_fa).exists():
                sel = ad / "selected_candidate_contigs.fasta"
                kept = []
                for s in SeqIO.parse(str(asm_fa), "fasta"):
                    if len(s.seq) >= min_len:
                        kept.append(s)
                if kept:
                    SeqIO.write(kept, str(sel), "fasta")
                    for k in kept:
                        if len(k.seq) > best_len:
                            best_len = len(k.seq); best_contig = k.id
                    # blastn
                    b6 = bd / "blastn_vs_refined_mitogenome.tsv"
                    run_cmd(["blastn","-query",str(sel),"-db",str(refined_fa),"-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen","-evalue",str(ca.get('blastn_evalue','1e-10')),"-num_threads",str(threads),"-out",str(b6)], check=False)
                    if b6.exists() and b6.stat().st_size>0:
                        bdf = pd.read_csv(b6, sep="\t", header=None)
                        x = bdf.iloc[0]; bident = float(x[2]); qlen = float(x[12]); bqcov = round(100*float(x[3])/max(qlen,1),2); bsstart = int(x[8]); bsend = int(x[9])
                    # tblastn
                    refaa = _reference_protein(config, gene, code)
                    if refaa:
                        qfaa = bd / f"reference_{gene}.faa"
                        with open(qfaa,"w",encoding="utf-8") as f: f.write(f">{gene}\n{refaa}\n")
                        run_cmd(["makeblastdb","-in",str(sel),"-dbtype","nucl"], check=False)
                        t6 = bd / "tblastn_reference_protein_vs_candidate_assembly.tsv"
                        run_cmd(["tblastn","-query",str(qfaa),"-db",str(sel),"-outfmt","6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen","-evalue",str(ca.get('tblastn_evalue','1e-5')),"-num_threads",str(threads),"-out",str(t6)], check=False)
                        if t6.exists() and t6.stat().st_size>0:
                            tdf = pd.read_csv(t6, sep="\t", header=None)
                            y=tdf.iloc[0]; tbpid=float(y[2]); qlen=float(y[12]); slen=float(y[13]); al=float(y[3]); tbqcov=round(100*al/max(qlen,1),2); tbscov=round(100*al/max(slen,1),2); tbbs=float(y[11])
                    # recommendation
                    if tbpid != "." and tbqcov != "." and tbpid >= float(ca.get('tblastn_min_identity',25)) and tbqcov >= float(ca.get('tblastn_min_coverage',25)):
                        recm = "LOCAL_ASSEMBLY_SUPPORTS_MISSING_GENE"
                        if bqcov != "." and bqcov < float(ca.get('blastn_min_coverage',30)):
                            recm = "LOCAL_ASSEMBLY_SUGGESTS_MITOGENOME_ERROR"
                    elif tbpid != ".":
                        recm = "LOCAL_ASSEMBLY_NEEDS_REVIEW"
            row = {"gene":gene,"target_id":tid,"read_set":rs,"assembler":assembler or ".","assembly_fasta":str(asm_fa) if asm_fa else ".","best_contig_id":best_contig,"best_contig_len":best_len,"blastn_best_hit_refined_seqid":rec.id if bident != "." else ".","blastn_best_pident":bident,"blastn_best_qcov":bqcov,"blastn_best_sstart":bsstart,"blastn_best_send":bsend,"tblastn_best_pident":tbpid,"tblastn_best_query_coverage":tbqcov,"tblastn_best_subject_coverage":tbscov,"tblastn_best_bitscore":tbbs,"reference_protein_length":len(_reference_protein(config,gene,code) or ""),"candidate_region_status":status,"recommendation":recm,"comment":"diagnostic-only"}
            rows.append(row)
            pd.DataFrame([row]).to_csv(dd / "candidate_gene_diagnosis.tsv", sep="\t", index=False)

    sdf = pd.DataFrame(rows)
    sdf.to_csv(sum_tsv, sep="\t", index=False)
    with open(sum_md, "w", encoding="utf-8") as md:
        md.write("# Candidate assembly summary\n\n")
        if sdf.empty:
            md.write("No candidate assemblies generated.\n")
        else:
            for g, grp in sdf.groupby("gene"):
                best = grp.iloc[0]
                md.write(f"## {g}\n\n")
                md.write(f"- best candidate: {best['target_id']}\n")
                md.write(f"- read_sets used: {','.join(sorted(set(grp['read_set'].astype(str))))}\n")
                md.write(f"- contigs montados: {int((grp['best_contig_len']>0).sum())}\n")
                md.write(f"- melhor tblastn pident: {best['tblastn_best_pident']}\n")
                md.write(f"- melhor blastn pident: {best['blastn_best_pident']}\n")
                md.write(f"- recomendação: {best['recommendation']}\n\n")
        md.write("Nenhum GenBank foi alterado nesta etapa.\n")
    return outdir
