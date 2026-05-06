from __future__ import annotations
from pathlib import Path
from .utils import ensure_dir, run_cmd, which_or_path


def _build_command(mf: dict, fasta_input: str, outdir: Path):
    mode = mf.get("mode", "python_interpreter")
    prefix = mf.get("prefix", "mitofinder")
    extra_args = mf.get("extra_args", []) or []

    # Conservative generic invocation that can be overridden by wrapper mode.
    if mode == "python_interpreter":
        py2 = mf.get("python2")
        script = mf.get("script")
        if not py2 or not Path(py2).exists():
            raise RuntimeError("MitoFinder mode=python_interpreter requires a valid tools.mitofinder.python2 path")
        if not script or not Path(script).exists():
            raise RuntimeError("MitoFinder mode=python_interpreter requires a valid tools.mitofinder.script path")
        cmd = [py2, script, "-j", prefix, "-a", fasta_input, "-o", str(outdir)] + list(extra_args)
    elif mode == "conda_env":
        conda = mf.get("conda_executable", "conda")
        env = mf.get("conda_env")
        script = mf.get("script")
        conda_path = which_or_path(conda)
        if not conda_path:
            raise RuntimeError("MitoFinder mode=conda_env requires a valid conda executable")
        if not env:
            raise RuntimeError("MitoFinder mode=conda_env requires tools.mitofinder.conda_env")
        if not script or not Path(script).exists():
            raise RuntimeError("MitoFinder mode=conda_env requires a valid tools.mitofinder.script path")
        cmd = [conda_path, "run", "-n", env, script, "-j", prefix, "-a", fasta_input, "-o", str(outdir)] + list(extra_args)
    elif mode == "wrapper":
        wrapper = mf.get("wrapper")
        if not wrapper or not Path(wrapper).exists():
            raise RuntimeError("MitoFinder mode=wrapper requires a valid tools.mitofinder.wrapper path")
        cmd = [wrapper, fasta_input, str(outdir)] + list(extra_args)
    else:
        raise RuntimeError(f"Unsupported MitoFinder mode: {mode}")

    return cmd


def run_mitofinder_for_fasta(config: dict, fasta_input: str, outdir: Path) -> Path:
    ensure_dir(outdir)
    tools = config.get("tools", {}) or {}
    mf = tools.get("mitofinder", {}) or {}
    if not mf.get("enabled", False):
        raise RuntimeError("Input is FASTA but tools.mitofinder.enabled is false")

    cmd = _build_command(mf, fasta_input, outdir)
    run_cmd(cmd, log_file=outdir / "mitofinder.log", check=True)

    # Expected handoff file (can be overridden in config if needed)
    annotated_gb = mf.get("annotated_genbank")
    if annotated_gb:
        gb = Path(annotated_gb)
    else:
        gb = outdir / "mitofinder_annotated.gb"

    if not gb.exists():
        raise RuntimeError(
            f"MitoFinder finished but annotated GenBank not found: {gb}. "
            "Set tools.mitofinder.annotated_genbank to the expected output file."
        )

    with open(outdir / "mitofinder_command.txt", "w", encoding="utf-8") as f:
        f.write(" ".join(cmd) + "\n")

    return gb
