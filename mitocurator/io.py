from __future__ import annotations
from pathlib import Path
from Bio import SeqIO

GENBANK_EXT = {".gb", ".gbk", ".genbank"}
FASTA_EXT = {".fa", ".fasta", ".fna"}

def infer_format(path: str | Path) -> str:
    suffix = Path(path).suffix.lower()
    if suffix in GENBANK_EXT:
        return "genbank"
    if suffix in FASTA_EXT:
        return "fasta"
    raise ValueError(f"Cannot infer format from extension: {path}")

def read_record(path: str | Path):
    fmt = infer_format(path)
    return next(SeqIO.parse(str(path), fmt)), fmt

def write_record(record, path: str | Path, fmt: str):
    SeqIO.write(record, str(path), fmt)

def feature_name(feat) -> str:
    for key in ("gene", "product", "locus_tag", "note"):
        if key in feat.qualifiers:
            return feat.qualifiers[key][0]
    return "."

def normalize_gene_name(name: str) -> str:
    if not name:
        return "."
    x = name.strip()
    repl = {
        "COI": "COX1", "COXI": "COX1", "CO1": "COX1",
        "COII": "COX2", "COXII": "COX2", "CO2": "COX2",
        "COIII": "COX3", "COXIII": "COX3", "CO3": "COX3",
        "CYTB": "CYTB", "COB": "CYTB",
        "ATPASE6": "ATP6", "ATPASE8": "ATP8",
    }
    u = x.upper().replace("-", "").replace("_", "")
    return repl.get(u, x)
