from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from .utils import load_config, ensure_dir, safe_get
from .check_tools import check_tools
from .gene_qc import diagnose
from .rotate import rotate_to_gene
from .refinement import refine_annotation
from .local_consensus import repair_cds_local_consensus
from .io import read_record, write_record
from .final_molecule_preparation import run_final_molecule_preparation
from .polish import run_global_polish


def outdir_from_config(config: dict) -> Path:
    outdir = safe_get(config, ["output", "outdir"], None)
    if not outdir:
        project = safe_get(config, ["project", "name"], "mitocurator_run")
        outdir = str(Path.cwd() / project)
    return ensure_dir(outdir)


def _load_problems(ref_dir: Path) -> list[dict]:
    """Build problems list from refinement TSVs for repair_cds_local_consensus.

    Sources (never overlap — cds_tsv covers annotated genes, miss_tsv covers absent genes):
    - cds_refinement_candidates.tsv: all rows (INTERNAL_STOP, INCOMPLETE_LENGTH,
      INCOMPLETE_COVERAGE). Column aliases old_start/old_end/old_strand/
      old_internal_stop_count are resolved inside repair_cds_local_consensus.
    - missing_gene_candidates.tsv: top-1 candidate per gene (TSV is already
      sorted by internal_stop_count ASC, length_nt DESC).
    """
    problems: list[dict] = []

    cds_tsv = ref_dir / "cds_refinement_candidates.tsv"
    if cds_tsv.exists():
        df = pd.read_csv(cds_tsv, sep="\t")
        if not df.empty:
            problems.extend(df.to_dict("records"))

    miss_tsv = ref_dir / "missing_gene_candidates.tsv"
    if miss_tsv.exists():
        df = pd.read_csv(miss_tsv, sep="\t")
        if not df.empty:
            problems.extend(df.drop_duplicates(subset="gene", keep="first").to_dict("records"))

    return problems


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


def cmd_polish(args):
    config = load_config(args.config)
    assembly = Path(args.input) if args.input else Path(
        safe_get(config, ["polish", "assembly_fasta"], "") or ""
    )
    if not assembly or not assembly.exists():
        raise SystemExit(
            "Assembly FASTA not found. Provide --input or set polish.assembly_fasta in config."
        )
    root = outdir_from_config(config)
    out = run_global_polish(config, assembly, root)
    print(f"Polished FASTA written to: {out}")


def cmd_final_molecule(args):
    config = load_config(args.config)
    input_gb = Path(config["input"]["mitogenome"])
    step_name = safe_get(config, ["output", "step_dirs", "final_molecule"], "08_final_molecule")
    outdir = ensure_dir(outdir_from_config(config) / step_name)
    out = run_final_molecule_preparation(config, input_gb, outdir)
    print(f"Final molecule written to: {out}")


def cmd_run(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    print(f"MitoCurator run directory: {root}")

    logs = ensure_dir(root / "00_logs")
    tc = check_tools(config, logs)
    print(f"[1/8] Tool check: {tc}")

    # [2/8] Global polish — optional; runs on raw assembly FASTA before annotation.
    # MitoFinder must be re-run on the polished FASTA (external step) before cmd_run
    # can proceed to refinement. If polish is disabled or assembly_fasta is absent,
    # this step is skipped and cmd_run continues from config["input"]["mitogenome"].
    polish_enabled = bool(safe_get(config, ["polish", "enabled"], False))
    polish_dir = None
    if polish_enabled:
        assembly_fasta_path = safe_get(config, ["polish", "assembly_fasta"], None)
        if assembly_fasta_path and Path(assembly_fasta_path).exists():
            polish_dir = root
            polished = run_global_polish(config, Path(assembly_fasta_path), root)
            print(f"[2/8] Polish: {polished}")
            print(
                "      Next: annotate the polished FASTA with MitoFinder, then set "
                "input.mitogenome to the resulting GenBank and re-run."
            )
        else:
            print("[2/8] Polish: skipped (polish.assembly_fasta not set or not found)")
    else:
        print("[2/8] Polish: disabled")

    annotated_gb = config["input"]["mitogenome"]
    print(f"[3/8] MitoFinder annotation: {annotated_gb}")

    refinement_enabled = bool(safe_get(config, ["refinement", "enabled"], True))
    refined_gb = annotated_gb
    ref_dir = None
    if refinement_enabled:
        ref_dir = ensure_dir(root / "05_refinement")
        refined_gb = refine_annotation(config, annotated_gb, ref_dir)
        print(f"[4/8] Annotation refinement: {refined_gb}")
    else:
        print("[4/8] Annotation refinement: disabled")

    # [5/8] Local consensus repair (suggest → review audit_log → re-run with apply).
    # mode=suggest: only candidate FASTAs and audit_log written; record unchanged.
    # mode=apply:   record modified in-place; repaired.gb written and passed downstream.
    lc_enabled = bool(safe_get(config, ["local_consensus", "enabled"], True))
    lc_mode = safe_get(config, ["local_consensus", "mode"], "suggest")
    lc_gb = refined_gb
    lc_dir = None
    if lc_enabled:
        if ref_dir is None:
            print("[5/8] Local consensus: skipped (refinement disabled — TSVs unavailable)")
        else:
            problems = _load_problems(ref_dir)
            if not problems:
                print("[5/8] Local consensus: skipped (no problems found in refinement TSVs)")
            else:
                step_name = safe_get(config, ["output", "step_dirs", "local_consensus"],
                                     "06_local_consensus")
                lc_dir = ensure_dir(root / step_name)
                record, fmt = read_record(refined_gb)
                repair_cds_local_consensus(config, record, problems, root, mode=lc_mode)
                if lc_mode == "apply":
                    lc_gb = lc_dir / "repaired.gb"
                    write_record(record, lc_gb, fmt)
                print(f"[5/8] Local consensus ({lc_mode}): {lc_dir}")
    else:
        print("[5/8] Local consensus: disabled")

    try:
        rot_dir = ensure_dir(root / "04_rotation")
        config["input"]["mitogenome"] = str(lc_gb)
        rotated_input = rotate_to_gene(config, rot_dir)
        print(f"[6/8] Rotation: {rotated_input}")
        config["input"]["mitogenome"] = str(rotated_input)
    except Exception as e:
        print(f"[6/8] Rotation skipped/failed: {e}")
        print("      Proceeding with current annotation for diagnosis.")

    qc_dir = ensure_dir(root / "07_gene_qc")
    diagnose(config, qc_dir)
    print(f"[7/8] Diagnosis: {qc_dir}")

    fm_dir = None
    fm_enabled = bool(safe_get(config, ["final_molecule", "enabled"], True))
    if fm_enabled:
        fm_step = safe_get(config, ["output", "step_dirs", "final_molecule"], "08_final_molecule")
        fm_dir = ensure_dir(root / fm_step)
        fm_input = Path(config["input"]["mitogenome"])
        fm_cfg = safe_get(config, ["final_molecule"], {}) or {}
        fm_out = run_final_molecule_preparation(
            config,
            fm_input,
            fm_dir,
            at_window_size=int(fm_cfg.get("at_window_size", 500)),
            at_min_pct=float(fm_cfg.get("at_min_pct", 85.0)),
            at_min_len=int(fm_cfg.get("at_min_len", 500)),
            expected_length=fm_cfg.get("expected_length"),
        )
        print(f"[8/8] Final molecule: {fm_out}")
    else:
        print("[8/8] Final molecule preparation: disabled")

    print("\nMain outputs:")
    print(f"  {logs / 'tool_check.tsv'}")
    if polish_dir is not None:
        step_name = safe_get(config, ["output", "step_dirs", "polish"], "01_polish")
        print(f"  {root / step_name / 'audit_log.jsonl'}")
        print(f"  {root / step_name / 'polished.fasta'}")
    print(f"  {root / '05_refinement' / 'refined.gb'}")
    print(f"  {root / '05_refinement' / 'expected_gene_set.tsv'}")
    print(f"  {root / '05_refinement' / 'added_features.tsv'}")
    print(f"  {root / '05_refinement' / 'missing_gene_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'cds_refinement_candidates.tsv'}")
    if lc_dir is not None:
        print(f"  {lc_dir / 'audit_log.jsonl'}")
        print(f"  {lc_dir / 'summary.tsv'}")
        if lc_mode == "apply":
            print(f"  {lc_dir / 'repaired.gb'}")
    print(f"  {qc_dir / 'gene_qc.tsv'}")
    print(f"  {qc_dir / 'problematic_features.tsv'}")
    print(f"  {qc_dir / 'intergenic_regions.tsv'}")
    print(f"  {qc_dir / 'diagnostic_summary.md'}")
    if fm_dir is not None:
        print(f"  {fm_dir / 'final_molecule.gb'}")
        print(f"  {fm_dir / 'final_molecule_report.md'}")

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

    p_pol = sub.add_parser("polish", help="Polish raw assembly FASTA with HiFi and/or Illumina reads")
    p_pol.add_argument("--config", required=True)
    p_pol.add_argument("--input", default=None,
                       help="Raw assembly FASTA (overrides polish.assembly_fasta in config)")
    p_pol.set_defaults(func=cmd_polish)

    p_fm = sub.add_parser("final-molecule", help="Prepare final GenBank/FASTA and annotate A+T-rich region")
    p_fm.add_argument("--config", required=True)
    p_fm.set_defaults(func=cmd_final_molecule)

    return p

def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
