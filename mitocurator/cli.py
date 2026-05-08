from __future__ import annotations
import argparse
from pathlib import Path
from .utils import load_config, ensure_dir, safe_get
from .check_tools import check_tools
from .gene_qc import diagnose
from .rotate import rotate_to_gene
from .refinement import refine_annotation
from .read_support import run_read_support
from .targeted_extraction import run_targeted_extraction
from .reconstruction_pools import run_reconstruction_pools
from .targeted_consensus import run_targeted_consensus

def outdir_from_config(config: dict) -> Path:
    outdir = safe_get(config, ["output", "outdir"], None)
    if not outdir:
        base = safe_get(config, ["project", "output_base_dir"], None)
        prefix = safe_get(config, ["project", "output_prefix"], None) or safe_get(config, ["project", "name"], "mitocurator_run")
        outdir = str(Path(base) / str(prefix)) if base else str(Path.cwd() / str(prefix))
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

def cmd_read_support(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    run_read_support(config, Path(safe_get(config, ["input", "mitogenome"])), root / "05_refinement", ensure_dir(root / "06_read_support"))
    print(f"Read-support files written to: {root / '06_read_support'}")

def cmd_targeted_extraction(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    run_targeted_extraction(config, root, root / "05_refinement", root / "06_read_support", ensure_dir(root / "08_targeted_extraction"))
    print(f"Targeted extraction files written to: {root / '08_targeted_extraction'}")

def cmd_reconstruction_pools(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    run_reconstruction_pools(config, root, root / "06_read_support", root / "08_targeted_extraction", ensure_dir(root / "09_reconstruction_pools"))
    print(f"Reconstruction pools written to: {root / '09_reconstruction_pools'}")

def cmd_targeted_consensus(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    run_targeted_consensus(config, root, root / "05_refinement", root / "09_reconstruction_pools", ensure_dir(root / "10_targeted_consensus"))
    print(f"Targeted consensus written to: {root / '10_targeted_consensus'}")

def cmd_run(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    print(f"MitoCurator run directory: {root}")

    logs = ensure_dir(root / "00_logs")
    tc = check_tools(config, logs)
    print(f"[1/9] Tool check: {tc}")

    annotated_gb = config["input"]["mitogenome"]
    print(f"[2/9] MitoFinder annotation: {annotated_gb}")

    refinement_enabled = bool(safe_get(config, ["refinement", "enabled"], True))
    refined_gb = annotated_gb
    if refinement_enabled:
        ref_dir = ensure_dir(root / "05_refinement")
        refined_gb = refine_annotation(config, annotated_gb, ref_dir)
        print(f"[3/9] Annotation refinement: {refined_gb}")
    else:
        print("[3/9] Annotation refinement: disabled")

    try:
        rot_dir = ensure_dir(root / "04_rotation")
        config["input"]["mitogenome"] = str(refined_gb)
        rotated_input = rotate_to_gene(config, rot_dir)
        print(f"[4/9] Rotation: {rotated_input}")
        config["input"]["mitogenome"] = str(rotated_input)
    except Exception as e:
        print(f"[4/9] Rotation skipped/failed: {e}")
        print("      Proceeding with current annotation for diagnosis.")

    if bool(safe_get(config, ["read_support", "enabled"], True)):
        rs_dir = ensure_dir(root / "06_read_support")
        run_read_support(config, Path(refined_gb), root / "05_refinement", rs_dir)
        print(f"[5/9] Read support: {rs_dir}")
    else:
        print("[5/9] Read support: disabled")

    if bool(safe_get(config, ["targeted_extraction", "enabled"], True)):
        te_dir = ensure_dir(root / "08_targeted_extraction")
        run_targeted_extraction(config, root, root / "05_refinement", root / "06_read_support", te_dir)
        print(f"[6/9] Targeted extraction: {te_dir}")
    else:
        print("[6/9] Targeted extraction: disabled")

    if bool(safe_get(config, ["reconstruction_pools", "enabled"], True)):
        rp_dir = ensure_dir(root / "09_reconstruction_pools")
        run_reconstruction_pools(config, root, root / "06_read_support", root / "08_targeted_extraction", rp_dir)
        print(f"[7/9] Reconstruction pools: {rp_dir}")
    else:
        print("[7/9] Reconstruction pools: disabled")

    if bool(safe_get(config, ["targeted_consensus", "enabled"], True)):
        tc_dir = ensure_dir(root / "10_targeted_consensus")
        run_targeted_consensus(config, root, root / "05_refinement", root / "09_reconstruction_pools", tc_dir)
        print(f"[8/9] Targeted consensus: {tc_dir}")
    else:
        print("[8/9] Targeted consensus: disabled")

    qc_dir = ensure_dir(root / "07_gene_qc")
    diagnose(config, qc_dir)
    print(f"[9/9] Diagnosis: {qc_dir}")

    print("\nMain outputs:")
    print(f"  {logs / 'tool_check.tsv'}")
    print(f"  {root / '05_refinement' / 'refined.gb'}")
    print(f"  {root / '05_refinement' / 'expected_gene_set.tsv'}")
    print(f"  {root / '05_refinement' / 'added_features.tsv'}")
    print(f"  {root / '05_refinement' / 'missing_gene_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'cds_refinement_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'reference_similarity_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_reference_check.tsv'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_proteins.faa'}")
    print(f"  {root / '05_refinement' / 'curation_recommendations.tsv'}")
    print(f"  {root / '05_refinement' / 'curation_recommendations.md'}")
    print(f"  {root / '05_refinement' / 'missing_gene_candidate_proteins.faa'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_reference_alignment.tsv'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_stop_context.tsv'}")
    print(f"  {qc_dir / 'gene_qc.tsv'}")
    print(f"  {qc_dir / 'problematic_features.tsv'}")
    print(f"  {qc_dir / 'intergenic_regions.tsv'}")
    print(f"  {qc_dir / 'diagnostic_summary.md'}")
    print(f"  {root / '06_read_support' / 'readset_consensus_recommendations.tsv'}")
    print(f"  {root / '06_read_support' / 'readset_consensus_recommendations.md'}")
    print(f"  {root / '08_targeted_extraction' / 'targeted_read_extraction.tsv'}")
    print(f"  {root / '08_targeted_extraction' / 'targeted_read_extraction.md'}")
    print(f"  {root / '08_targeted_extraction' / 'targets.bed'}")
    print(f"  {root / '09_reconstruction_pools' / 'reconstruction_pools.tsv'}")
    print(f"  {root / '09_reconstruction_pools' / 'reconstruction_pools.md'}")
    print(f"  {root / '10_targeted_consensus' / 'targeted_consensus.tsv'}")
    print(f"  {root / '10_targeted_consensus' / 'targeted_consensus.md'}")

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

    for name, help_txt, func in [
        ("read-support", "Run read support stage only", cmd_read_support),
        ("targeted-extraction", "Run targeted extraction stage only", cmd_targeted_extraction),
        ("reconstruction-pools", "Run reconstruction pools stage only", cmd_reconstruction_pools),
        ("targeted-consensus", "Run targeted consensus stage only", cmd_targeted_consensus),
    ]:
        p = sub.add_parser(name, help=help_txt)
        p.add_argument("--config", required=True)
        p.set_defaults(func=func)

    return p

def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
