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
from .read_support import run_read_support
from .targeted_extraction import run_targeted_extraction
from .reconstruction_pools import run_reconstruction_pools
from .targeted_consensus import run_targeted_consensus
from .candidate_assembly import run_candidate_assembly


def outdir_from_config(config: dict) -> Path:
    legacy_outdir = safe_get(config, ["output", "outdir"], None)
    if legacy_outdir:
        return ensure_dir(legacy_outdir)

    output_base_dir = safe_get(config, ["project", "output_base_dir"], None)
    output_prefix = safe_get(config, ["project", "output_prefix"], None)

    if output_base_dir and output_prefix:
        return ensure_dir(Path(output_base_dir) / output_prefix)

    project = safe_get(config, ["project", "name"], "mitocurator_run")
    return ensure_dir(Path.cwd() / project)


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
    refined_gb = root / "05_refinement" / "refined.gb"
    if not refined_gb.exists():
        raise FileNotFoundError(f"Missing refined GenBank: {refined_gb}")
    rs_dir = ensure_dir(root / "06_read_support")
    run_read_support(config, refined_gb, root / "05_refinement", rs_dir)
    print(f"Read support written to: {rs_dir}")


def cmd_targeted_extraction(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    te_dir = ensure_dir(root / "08_targeted_extraction")
    run_targeted_extraction(
        config,
        root,
        root / "05_refinement",
        root / "06_read_support",
        te_dir,
    )
    print(f"Targeted extraction written to: {te_dir}")


def cmd_reconstruction_pools(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    pools_dir = ensure_dir(root / "09_reconstruction_pools")
    run_reconstruction_pools(
        config,
        root,
        root / "06_read_support",
        root / "08_targeted_extraction",
        pools_dir,
    )
    print(f"Reconstruction pools written to: {pools_dir}")


def cmd_targeted_consensus(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    cons_dir = ensure_dir(root / "10_targeted_consensus")
    run_targeted_consensus(
        config,
        root,
        root / "05_refinement",
        root / "09_reconstruction_pools",
        cons_dir,
    )
    print(f"Targeted consensus written to: {cons_dir}")


def cmd_candidate_assembly(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    ca_dir = ensure_dir(root / "11_candidate_assembly")
    run_candidate_assembly(
        config,
        root,
        root / "10_targeted_consensus",
        root / "09_reconstruction_pools",
        root / "05_refinement",
        ca_dir,
    )
    print(f"Candidate assembly written to: {ca_dir}")


def cmd_run(args):
    config = load_config(args.config)
    root = outdir_from_config(config)
    print(f"MitoCurator run directory: {root}")

    logs = ensure_dir(root / "00_logs")
    tc = check_tools(config, logs)
    print(f"[1/10] Tool check: {tc}")

    current_input = config["input"]["mitogenome"]
    fmt = infer_format(current_input)

    if fmt == "fasta":
        mf_dir = ensure_dir(root / "03_mitofinder")
        annotated_gb = run_mitofinder_for_fasta(config, current_input, mf_dir)
        config["input"]["mitogenome"] = str(annotated_gb)
        print(f"[2/10] MitoFinder annotation: {annotated_gb}")
    else:
        annotated_gb = current_input
        config["input"]["mitogenome"] = str(annotated_gb)
        print("[2/10] MitoFinder annotation: skipped (input already annotated GenBank)")

    refinement_enabled = bool(safe_get(config, ["refinement", "enabled"], True))
    refined_gb = annotated_gb

    if refinement_enabled:
        ref_dir = ensure_dir(root / "05_refinement")
        refined_gb = refine_annotation(config, annotated_gb, ref_dir)
        config["input"]["mitogenome"] = str(refined_gb)
        print(f"[3/10] Annotation refinement: {refined_gb}")
    else:
        print("[3/10] Annotation refinement: disabled")

    try:
        rot_dir = ensure_dir(root / "04_rotation")
        config["input"]["mitogenome"] = str(refined_gb)
        rotated_input = rotate_to_gene(config, rot_dir)
        config["input"]["mitogenome"] = str(rotated_input)
        print(f"[4/10] Rotation: {rotated_input}")
    except Exception as e:
        print(f"[4/10] Rotation skipped/failed: {e}")
        print("      Proceeding with current annotation for downstream steps.")
        config["input"]["mitogenome"] = str(refined_gb)

    read_support_enabled = bool(safe_get(config, ["read_support", "enabled"], False))
    if read_support_enabled:
        rs_dir = ensure_dir(root / "06_read_support")
        run_read_support(config, Path(refined_gb), root / "05_refinement", rs_dir)
        print(f"[5/10] Read support: {rs_dir}")
    else:
        print("[5/10] Read support: disabled")

    targeted_enabled = bool(safe_get(config, ["targeted_extraction", "enabled"], False))
    if targeted_enabled:
        te_dir = ensure_dir(root / "08_targeted_extraction")
        run_targeted_extraction(
            config,
            root,
            root / "05_refinement",
            root / "06_read_support",
            te_dir,
        )
        print(f"[6/10] Targeted extraction: {te_dir}")
    else:
        print("[6/10] Targeted extraction: disabled")

    pools_enabled = bool(safe_get(config, ["reconstruction_pools", "enabled"], False))
    if pools_enabled:
        pools_dir = ensure_dir(root / "09_reconstruction_pools")
        run_reconstruction_pools(
            config,
            root,
            root / "06_read_support",
            root / "08_targeted_extraction",
            pools_dir,
        )
        print(f"[7/10] Reconstruction pools: {pools_dir}")
    else:
        print("[7/10] Reconstruction pools: disabled")

    consensus_enabled = bool(safe_get(config, ["targeted_consensus", "enabled"], False))
    if consensus_enabled:
        cons_dir = ensure_dir(root / "10_targeted_consensus")
        run_targeted_consensus(
            config,
            root,
            root / "05_refinement",
            root / "09_reconstruction_pools",
            cons_dir,
        )
        print(f"[8/10] Targeted consensus: {cons_dir}")
    else:
        print("[8/10] Targeted consensus: disabled")

    candidate_assembly_enabled = bool(safe_get(config, ["candidate_assembly", "enabled"], False))
    if candidate_assembly_enabled:
        ca_dir = ensure_dir(root / "11_candidate_assembly")
        run_candidate_assembly(
            config,
            root,
            root / "10_targeted_consensus",
            root / "09_reconstruction_pools",
            root / "05_refinement",
            ca_dir,
        )
        print(f"[9/10] Candidate assembly: {ca_dir}")
    else:
        print("[9/10] Candidate assembly: disabled")

    qc_dir = ensure_dir(root / "07_gene_qc")
    diagnose(config, qc_dir)
    print(f"[10/10] Diagnosis: {qc_dir}")

    print("\nMain outputs:")
    print(f"  {logs / 'tool_check.tsv'}")
    print(f"  {root / '05_refinement' / 'refined.gb'}")
    print(f"  {root / '05_refinement' / 'expected_gene_set.tsv'}")
    print(f"  {root / '05_refinement' / 'added_features.tsv'}")
    print(f"  {root / '05_refinement' / 'missing_gene_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'cds_refinement_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'reference_similarity_candidates.tsv'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_reference_check.tsv'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_stop_context.tsv'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_reference_alignment.tsv'}")
    print(f"  {root / '05_refinement' / 'missing_gene_candidate_proteins.faa'}")
    print(f"  {root / '05_refinement' / 'problematic_cds_proteins.faa'}")
    print(f"  {root / '05_refinement' / 'curation_recommendations.tsv'}")
    print(f"  {root / '05_refinement' / 'curation_recommendations.md'}")

    if read_support_enabled:
        print(f"  {root / '06_read_support' / 'problematic_stop_read_support.tsv'}")
        print(f"  {root / '06_read_support' / 'problematic_stop_variants.tsv'}")
        print(f"  {root / '06_read_support' / 'read_support_summary.md'}")
        print(f"  {root / '06_read_support' / 'readset_consensus_recommendations.tsv'}")
        print(f"  {root / '06_read_support' / 'readset_consensus_recommendations.md'}")

    if targeted_enabled:
        print(f"  {root / '08_targeted_extraction' / 'targets.bed'}")
        print(f"  {root / '08_targeted_extraction' / 'targeted_read_extraction.tsv'}")
        print(f"  {root / '08_targeted_extraction' / 'targeted_read_extraction.md'}")

    if pools_enabled:
        print(f"  {root / '09_reconstruction_pools' / 'reconstruction_pools.tsv'}")
        print(f"  {root / '09_reconstruction_pools' / 'reconstruction_pools.md'}")

    if consensus_enabled:
        print(f"  {root / '10_targeted_consensus' / 'targeted_consensus.tsv'}")
        print(f"  {root / '10_targeted_consensus' / 'targeted_consensus.md'}")
        print(f"  {root / '10_targeted_consensus' / 'best_missing_gene_candidates.tsv'}")
        print(f"  {root / '10_targeted_consensus' / 'best_missing_gene_candidates.md'}")
        print(f"  {root / '10_targeted_consensus' / 'cross_readset_missing_gene_candidates.tsv'}")
        print(f"  {root / '10_targeted_consensus' / 'cross_readset_missing_gene_candidates.md'}")

    if candidate_assembly_enabled:
        print(f"  {root / '11_candidate_assembly' / 'candidate_assembly_targets.tsv'}")
        print(f"  {root / '11_candidate_assembly' / 'candidate_assembly_summary.tsv'}")
        print(f"  {root / '11_candidate_assembly' / 'candidate_assembly_summary.md'}")

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

    p_rs = sub.add_parser("read-support", help="Run only the read-support stage")
    p_rs.add_argument("--config", required=True)
    p_rs.set_defaults(func=cmd_read_support)

    p_te = sub.add_parser("targeted-extraction", help="Run only the targeted-extraction stage")
    p_te.add_argument("--config", required=True)
    p_te.set_defaults(func=cmd_targeted_extraction)

    p_rp = sub.add_parser("reconstruction-pools", help="Run only the reconstruction-pools stage")
    p_rp.add_argument("--config", required=True)
    p_rp.set_defaults(func=cmd_reconstruction_pools)

    p_tc = sub.add_parser("targeted-consensus", help="Run only the targeted-consensus stage")
    p_tc.add_argument("--config", required=True)
    p_tc.set_defaults(func=cmd_targeted_consensus)

    p_ca = sub.add_parser("candidate-assembly", help="Run only the candidate-assembly stage")
    p_ca.add_argument("--config", required=True)
    p_ca.set_defaults(func=cmd_candidate_assembly)

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
