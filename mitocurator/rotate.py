from __future__ import annotations
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from .io import read_record, write_record, feature_name, normalize_gene_name
from .utils import ensure_dir

def find_gene_feature(record, gene_name: str):
    target = normalize_gene_name(gene_name).upper()
    for feat in record.features:
        if feat.type == "source":
            continue
        names = []
        for key in ("gene", "product", "note"):
            names.extend(feat.qualifiers.get(key, []))
        normed = [normalize_gene_name(x).upper() for x in names]
        if target in normed:
            return feat
    return None

def rotate_seq(seq, cut_1based: int):
    # cut_1based is first base of new sequence
    i = cut_1based - 1
    return seq[i:] + seq[:i]

def remap_pos(pos0: int, cut0: int, n: int) -> int:
    # input/output 0-based
    return (pos0 - cut0) % n

def remap_feature_simple(feat, cut0: int, n: int):
    start0 = int(feat.location.start)
    end0 = int(feat.location.end)
    length = end0 - start0
    new_start = remap_pos(start0, cut0, n)
    new_end = new_start + length
    new_feat = feat
    if new_end <= n:
        new_feat.location = FeatureLocation(new_start, new_end, strand=feat.location.strand)
    else:
        # Feature crosses new origin; represent as join.
        part1 = FeatureLocation(new_start, n, strand=feat.location.strand)
        part2 = FeatureLocation(0, new_end - n, strand=feat.location.strand)
        new_feat.location = CompoundLocation([part1, part2])
    return new_feat

def rotate_to_gene(config: dict, outdir: Path):
    ensure_dir(outdir)
    mitogenome = config["input"]["mitogenome"]
    gene = config["project"]["rotate_to"]
    rec, fmt = read_record(mitogenome)
    feat = find_gene_feature(rec, gene)
    if feat is None:
        raise RuntimeError(f"Gene `{gene}` was not found in current annotation. In v0.1, rotation requires an annotated GenBank input.")
    cut_1based = int(feat.location.start) + 1
    cut0 = cut_1based - 1
    n = len(rec.seq)

    new_seq = rotate_seq(rec.seq, cut_1based)
    new = rec[:]
    new.seq = new_seq
    new.id = rec.id + f"_rotated_to_{gene}"
    new.name = new.id[:16]
    new.description = rec.description + f" [rotated to {gene}]"

    remapped = []
    for f in rec.features:
        if f.type == "source":
            nf = f
            nf.location = FeatureLocation(0, n, strand=1)
        else:
            nf = remap_feature_simple(f, cut0, n)
        remapped.append(nf)
    new.features = remapped

    out_gb = outdir / f"{Path(mitogenome).stem}.rotated_to_{gene}.gb"
    write_record(new, out_gb, "genbank")

    with open(outdir / "rotation_report.tsv", "w", encoding="utf-8") as r:
        r.write("rotate_to\told_start\told_end\tstrand\tcut_position_1based\toutput\n")
        strand = "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."
        r.write(f"{gene}\t{int(feat.location.start)+1}\t{int(feat.location.end)}\t{strand}\t{cut_1based}\t{out_gb}\n")
    return out_gb
