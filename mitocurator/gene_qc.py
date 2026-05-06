from __future__ import annotations
from pathlib import Path
from collections import defaultdict
from Bio.Seq import Seq
import pandas as pd
from .io import read_record, feature_name, normalize_gene_name
from .utils import ensure_dir

def feature_bounds(feat):
    return int(feat.location.start) + 1, int(feat.location.end)

def strand_symbol(feat):
    return "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."

def extract_feature_seq(record, feat):
    return feat.extract(record.seq)

def translate_cds(seq, genetic_code: int):
    # Do not trim terminal stop; user needs explicit reporting.
    return str(Seq(str(seq)).translate(table=genetic_code, to_stop=False))

def stop_report(aa: str, nt: str):
    stops = []
    for i, a in enumerate(aa, start=1):
        if a == "*":
            status = "terminal" if i == len(aa) else "internal"
            codon = nt[(i-1)*3:i*3]
            stops.append((i, codon, status))
    internal = [s for s in stops if s[2] == "internal"]
    terminal = any(s[2] == "terminal" for s in stops)
    return stops, internal, terminal

def summarize_features(record, genetic_code: int):
    rows = []
    for feat in record.features:
        if feat.type == "source":
            continue
        start, end = feature_bounds(feat)
        name = feature_name(feat)
        product = feat.qualifiers.get("product", ["."])[0]
        gene = feat.qualifiers.get("gene", [name])[0]
        norm = normalize_gene_name(gene if gene != "." else product)
        seq = extract_feature_seq(record, feat)
        length_nt = len(seq)
        row = {
            "seqid": record.id,
            "feature_type": feat.type,
            "gene": gene,
            "gene_normalized": norm,
            "product": product,
            "start": start,
            "end": end,
            "strand": strand_symbol(feat),
            "length_nt": length_nt,
            "length_aa": ".",
            "multiple_of_three": ".",
            "internal_stop_count": ".",
            "internal_stop_positions": ".",
            "terminal_stop": ".",
            "status": "OK",
            "decision_hint": "OK",
            "comment": ".",
        }
        if feat.type == "CDS":
            nt = str(seq).upper()
            aa = translate_cds(nt, genetic_code)
            stops, internal, terminal = stop_report(aa, nt)
            row["length_aa"] = len(aa)
            row["multiple_of_three"] = "yes" if length_nt % 3 == 0 else "no"
            row["internal_stop_count"] = len(internal)
            row["internal_stop_positions"] = ",".join([str(s[0]) for s in internal]) if internal else "."
            row["terminal_stop"] = "yes" if terminal else "no"
            if length_nt % 3 != 0 and internal:
                row["status"] = "PROBLEM"
                row["decision_hint"] = "CHECK_FRAMESHIFT_AND_INTERNAL_STOP"
            elif length_nt % 3 != 0:
                row["status"] = "PROBLEM"
                row["decision_hint"] = "CHECK_LENGTH_NOT_MULTIPLE_OF_THREE"
            elif internal:
                row["status"] = "PROBLEM"
                row["decision_hint"] = "CHECK_INTERNAL_STOP"
        rows.append(row)
    return pd.DataFrame(rows)

def intergenic_regions(record):
    intervals = []
    for feat in record.features:
        if feat.type == "source":
            continue
        start, end = feature_bounds(feat)
        intervals.append((start, end, feat.type, feature_name(feat)))
    intervals.sort()
    rows = []
    for i in range(len(intervals) - 1):
        a = intervals[i]
        b = intervals[i+1]
        gap_start = a[1] + 1
        gap_end = b[0] - 1
        if gap_end >= gap_start:
            seq = record.seq[gap_start-1:gap_end]
            length = len(seq)
            at = 100.0 * (str(seq).upper().count("A") + str(seq).upper().count("T")) / length if length else 0
            rows.append({
                "start": gap_start,
                "end": gap_end,
                "length": length,
                "AT_percent": round(at, 2),
                "upstream_feature": f"{a[2]}:{a[3]}",
                "downstream_feature": f"{b[2]}:{b[3]}",
            })
    return pd.DataFrame(rows)

def diagnose(config: dict, outdir: Path):
    ensure_dir(outdir)
    genetic_code = int(config["project"]["genetic_code"])
    mitogenome = config["input"]["mitogenome"]
    rec, fmt = read_record(mitogenome)

    feature_df = summarize_features(rec, genetic_code)
    intergenic_df = intergenic_regions(rec)

    feature_df.to_csv(outdir / "gene_qc.tsv", sep="\t", index=False)
    intergenic_df.to_csv(outdir / "intergenic_regions.tsv", sep="\t", index=False)

    problems = feature_df[feature_df["status"] != "OK"].copy()
    problems.to_csv(outdir / "problematic_features.tsv", sep="\t", index=False)

    # Simple summary markdown
    counts = feature_df["feature_type"].value_counts().to_dict() if not feature_df.empty else {}
    with open(outdir / "diagnostic_summary.md", "w", encoding="utf-8") as f:
        f.write("# MitoCurator diagnostic summary\n\n")
        f.write(f"- Input: `{mitogenome}`\n")
        f.write(f"- Input format: `{fmt}`\n")
        f.write(f"- Sequence ID: `{rec.id}`\n")
        f.write(f"- Length: `{len(rec.seq)}` bp\n")
        f.write(f"- Genetic code: `{genetic_code}`\n\n")
        f.write("## Feature counts\n\n")
        for k, v in sorted(counts.items()):
            f.write(f"- {k}: {v}\n")
        f.write("\n## Problematic features\n\n")
        if problems.empty:
            f.write("No problematic CDS detected by the current rules.\n")
        else:
            f.write(f"{len(problems)} problematic feature(s) detected. See `problematic_features.tsv`.\n")
        if not intergenic_df.empty:
            top = intergenic_df.sort_values("length", ascending=False).iloc[0]
            f.write("\n## Largest intergenic region\n\n")
            f.write(f"- Coordinates: {top['start']}..{top['end']}\n")
            f.write(f"- Length: {top['length']} bp\n")
            f.write(f"- A+T: {top['AT_percent']}%\n")
            f.write(f"- Between: {top['upstream_feature']} and {top['downstream_feature']}\n")
    return outdir
