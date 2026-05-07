from __future__ import annotations
import argparse
from pathlib import Path
from .utils import load_config, ensure_dir, safe_get
from .check_tools import check_tools
from .gene_qc import diagnose
from .rotate import rotate_to_gene
from .refinement import refine_annotation
from .io import infer_format
from .mitofinder_runner import run_mitofinder_for_fasta


def outdir_from_config(config: dict) -> Path:
    outdir = safe_get(config, ["output", "outdir"], None)
    if not outdir:
        project = safe_get(config, ["project", "name"], "mitocurator_run")
        outdir = str(Path.cwd() / project)
    return ensure_dir(outdir)


def cmd_check_tools(args):
    config = load_config(args.config)
    outdir = ensure_dir(outdir_from_config(config) / "00_logs")
    outfile = check_tools(config, outdir)
    print(f"Tool check written to: {outfile}")


def cmd_diagnose(args):
    config = load_config(args.config)
    outdir = ensure_dir(outdir_from_config(config) / "07_gene_qc")
    diagnose(config, outdir)
    print(f"Diagnostic files written to: {outdir}")


def cmd_rotate(args):
    config = load_config(args.config)
    outdir = ensure_dir(outdir_from_config(config) / "04_rotation")
    out = rotate_to_gene(config, outdir)
    print(f"Rotated GenBank written to: {out}")


def cmd_run(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    print(f"MitoCurator run directory: {root}")

    logs = ensure_dir(root / "00_logs")
    tc = check_tools(config, logs)
    print(f"[1/5] Tool check: {tc}")

    # FASTA-first flow: annotate with MitoFinder when the input is FASTA.
    current_input = config["input"]["mitogenome"]
    fmt = infer_format(current_input)

    if fmt == "fasta":
        mf_dir = ensure_dir(root / "03_mitofinder")
        annotated_gb = run_mitofinder_for_fasta(config, current_input, mf_dir)
        config["input"]["mitogenome"] = str(annotated_gb)
        print(f"[2/5] MitoFinder annotation: {annotated_gb}")
    else:
        annotated_gb = current_input
        config["input"]["mitogenome"] = str(annotated_gb)
        print("[2/5] MitoFinder annotation: skipped (input already annotated GenBank)")

    refinement_enabled = bool(safe_get(config, ["refinement", "enabled"], True))
    refined_gb = annotated_gb

    if refinement_enabled:
        ref_dir = ensure_dir(root / "05_refinement")
        refined_gb = refine_annotation(config, annotated_gb, ref_dir)
        config["input"]["mitogenome"] = str(refined_gb)
        print(f"[3/5] Annotation refinement: {refined_gb}")
    else:
        print("[3/5] Annotation refinement: disabled")

    try:
        rot_dir = ensure_dir(root / "04_rotation")
        config["input"]["mitogenome"] = str(refined_gb)
        rotated_input = rotate_to_gene(config, rot_dir)
        config["input"]["mitogenome"] = str(rotated_input)
        print(f"[4/5] Rotation: {rotated_input}")
    except Exception as e:
        print(f"[4/5] Rotation skipped/failed: {e}")
        print("      Proceeding with current annotation for diagnosis.")
        config["input"]["mitogenome"] = str(refined_gb)

    qc_dir = ensure_dir(root / "07_gene_qc")
    diagnose(config, qc_dir)
    print(f"[5/5] Diagnosis: {qc_dir}")

    print("\nMain outputs:")
    print(f"  {logs / 'tool_check.tsv'}")
    print(f"  {root / '05_refinement' / 'refined.gb'}")
    print(f"  {root / '05_refinement' / 'expected_gene_set.tsv'}")
    print(f"  {root / '05_refinement' / 'added_features.tsv'}")
    print(f"  {root / '05_refinement' / 'missing_gene_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'cds_refinement_candidates.tsv'}")
    print(f"  {qc_dir / 'gene_qc.tsv'}")
    print(f"  {qc_dir / 'problematic_features.tsv'}")
    print(f"  {qc_dir / 'intergenic_regions.tsv'}")
    print(f"  {qc_dir / 'diagnostic_summary.md'}")


def build_parser():
    p = argparse.ArgumentParser(prog="mitocurator", description="MitoCurator v0.1-dev")
    sub = p.add_subparsers(dest="command", required=True)

    p_check = sub.add_parser("check-tools", help="Check external tools configured in config.yaml")
    p_check.add_argument("--config", required=True)
    p_check.set_defaults(func=cmd_check_tools)

    p_diag = sub.add_parser("diagnose", help="Run gene-level diagnostic report")
    p_diag.add_argument("--config", required=True)
    p_diag.set_defaults(func=cmd_diagnose)

    p_rot = sub.add_parser("rotate", help="Rotate annotated GenBank to user-defined gene")
    p_rot.add_argument("--config", required=True)
    p_rot.set_defaults(func=cmd_rotate)

    p_run = sub.add_parser("run", help="Run initial all-in-one diagnostic workflow")
    p_run.add_argument("--config", required=True)
    p_run.set_defaults(func=cmd_run)

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
