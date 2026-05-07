from __future__ import annotations
from pathlib import Path
from Bio.Seq import Seq
import pandas as pd
from .io import read_record, feature_name, normalize_gene_name
from .utils import ensure_dir


FEATURE_COLUMNS = [
    "seqid",
    "feature_type",
    "gene",
    "gene_normalized",
    "product",
    "start",
    "end",
    "strand",
    "length_nt",
    "length_aa",
    "multiple_of_three",
    "internal_stop_count",
    "internal_stop_positions",
    "terminal_stop",
    "status",
    "decision_hint",
    "comment",
]


INTERGENIC_COLUMNS = [
    "start",
    "end",
    "length",
    "AT_percent",
    "upstream_feature",
    "downstream_feature",
]


def feature_bounds(feat):
    return int(feat.location.start) + 1, int(feat.location.end)


def strand_symbol(feat):
    return "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."


def extract_feature_seq(record, feat):
    return feat.extract(record.seq)


def translate_cds(seq, genetic_code: int):
    return str(Seq(str(seq)).translate(table=genetic_code, to_stop=False))


def stop_report(aa: str, nt: str):
    stops = []
    for i, a in enumerate(aa, start=1):
        if a == "*":
            status = "terminal" if i == len(aa) else "internal"
            codon = nt[(i - 1) * 3:i * 3]
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

    return pd.DataFrame(rows, columns=FEATURE_COLUMNS)


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
        b = intervals[i + 1]

        gap_start = a[1] + 1
        gap_end = b[0] - 1

        if gap_end >= gap_start:
            seq = record.seq[gap_start - 1:gap_end]
            length = len(seq)
            at = 100.0 * (
                str(seq).upper().count("A") + str(seq).upper().count("T")
            ) / length if length else 0

            rows.append({
                "start": gap_start,
                "end": gap_end,
                "length": length,
                "AT_percent": round(at, 2),
                "upstream_feature": f"{a[2]}:{a[3]}",
                "downstream_feature": f"{b[2]}:{b[3]}",
            })

    return pd.DataFrame(rows, columns=INTERGENIC_COLUMNS)


def basic_sequence_summary(record):
    seq = str(record.seq).upper()
    n = len(seq)
    gc = 100.0 * (seq.count("G") + seq.count("C")) / n if n else 0
    at = 100.0 * (seq.count("A") + seq.count("T")) / n if n else 0

    return {
        "seqid": record.id,
        "length_bp": n,
        "GC_percent": round(gc, 2),
        "AT_percent": round(at, 2),
        "N_count": seq.count("N"),
    }


def diagnose(config: dict, outdir: Path):
    ensure_dir(outdir)

    genetic_code = int(config["project"]["genetic_code"])
    mitogenome = config["input"]["mitogenome"]

    rec, fmt = read_record(mitogenome)

    summary = basic_sequence_summary(rec)
    feature_df = summarize_features(rec, genetic_code)
    intergenic_df = intergenic_regions(rec)

    feature_df.to_csv(outdir / "gene_qc.tsv", sep="\t", index=False)
    intergenic_df.to_csv(outdir / "intergenic_regions.tsv", sep="\t", index=False)

    if feature_df.empty:
        problems = pd.DataFrame(columns=FEATURE_COLUMNS)
    else:
        problems = feature_df[feature_df["status"] != "OK"].copy()

    problems.to_csv(outdir / "problematic_features.tsv", sep="\t", index=False)

    with open(outdir / "sequence_summary.tsv", "w", encoding="utf-8") as f:
        f.write("seqid\tlength_bp\tGC_percent\tAT_percent\tN_count\n")
        f.write(
            f"{summary['seqid']}\t{summary['length_bp']}\t"
            f"{summary['GC_percent']}\t{summary['AT_percent']}\t{summary['N_count']}\n"
        )

    counts = feature_df["feature_type"].value_counts().to_dict() if not feature_df.empty else {}

    with open(outdir / "diagnostic_summary.md", "w", encoding="utf-8") as f:
        f.write("# MitoCurator diagnostic summary\n\n")
        f.write(f"- Input: `{mitogenome}`\n")
        f.write(f"- Input format: `{fmt}`\n")
        f.write(f"- Sequence ID: `{rec.id}`\n")
        f.write(f"- Length: `{len(rec.seq)}` bp\n")
        f.write(f"- GC content: `{summary['GC_percent']}`%\n")
        f.write(f"- AT content: `{summary['AT_percent']}`%\n")
        f.write(f"- N count: `{summary['N_count']}`\n")
        f.write(f"- Genetic code: `{genetic_code}`\n\n")

        f.write("## Feature counts\n\n")
        if counts:
            for k, v in sorted(counts.items()):
                f.write(f"- {k}: {v}\n")
        else:
            f.write(
                "No annotated features were found. "
                "This usually means that the input is FASTA or an unannotated GenBank file.\n"
            )

        f.write("\n## Problematic features\n\n")
        if feature_df.empty:
            f.write(
                "No feature-level diagnosis was possible because no annotated features were found.\n"
            )
        elif problems.empty:
            f.write("No problematic CDS detected by the current rules.\n")
        else:
            f.write(f"{len(problems)} problematic feature(s) detected. See `problematic_features.tsv`.\n")

        f.write("\n## Intergenic regions\n\n")
        if intergenic_df.empty:
            f.write(
                "No intergenic regions could be calculated because fewer than two annotated features were found.\n"
            )
        else:
            top = intergenic_df.sort_values("length", ascending=False).iloc[0]
            f.write(f"- Largest interval: {top['start']}..{top['end']}\n")
            f.write(f"- Length: {top['length']} bp\n")
            f.write(f"- A+T: {top['AT_percent']}%\n")
            f.write(f"- Between: {top['upstream_feature']} and {top['downstream_feature']}\n")

    return outdir
