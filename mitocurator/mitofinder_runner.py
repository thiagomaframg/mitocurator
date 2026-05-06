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

def run_mitofinder_for_fasta(config, input_fasta, outdir):
    from pathlib import Path
    import shutil

    from .utils import run_cmd

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mf_cfg = config.get("mitofinder", {})

    if not mf_cfg.get("enabled", True):
        raise RuntimeError("MitoFinder is disabled in the config file.")

    mode = mf_cfg.get("mode", "python_interpreter")

    reference_gb = mf_cfg.get("reference_gb")
    organism_code = str(mf_cfg.get("organism_code", 5))
    threads = str(mf_cfg.get("threads", 8))
    memory_gb = str(mf_cfg.get("memory_gb", 40))

    if not reference_gb:
        raise ValueError("Missing mitofinder.reference_gb in config file.")

    run_dir = outdir / "run"
    run_dir.mkdir(parents=True, exist_ok=True)

    if mode == "python_interpreter":
        python2 = mf_cfg.get("python2", "/usr/bin/python2.7")
        script = mf_cfg.get("script", "/home/thiagomafra/bin/MitoFinder/mitofinder")

        cmd = [
            python2,
            script,
            "-j", "mitofinder",
            "-a", str(input_fasta),
            "-r", str(reference_gb),
            "-o", organism_code,
            "-p", threads,
            "-m", memory_gb,
            "--override",
        ]

    elif mode == "wrapper":
        wrapper = mf_cfg.get("wrapper")
        if not wrapper:
            raise ValueError("Missing mitofinder.wrapper in config file.")

        cmd = [
            wrapper,
            str(input_fasta),
            str(reference_gb),
            organism_code,
            threads,
            memory_gb,
        ]

    elif mode == "conda_env":
        conda_executable = mf_cfg.get("conda_executable")
        conda_env = mf_cfg.get("conda_env")
        script = mf_cfg.get("script", "/home/thiagomafra/bin/MitoFinder/mitofinder")

        if not conda_executable:
            raise ValueError("Missing mitofinder.conda_executable in config file.")
        if not conda_env:
            raise ValueError("Missing mitofinder.conda_env in config file.")

        cmd = [
            conda_executable,
            "run",
            "-n", conda_env,
            "python2",
            script,
            "-j", "mitofinder",
            "-a", str(input_fasta),
            "-r", str(reference_gb),
            "-o", organism_code,
            "-p", threads,
            "-m", memory_gb,
            "--override",
        ]

    else:
        raise ValueError(f"Unsupported MitoFinder mode: {mode}")

    run_cmd(cmd, cwd=run_dir, log_file=outdir / "mitofinder.log", check=True)

    final_dir = run_dir / "mitofinder" / "mitofinder_MitoFinder_mitfi_Final_Results"

    if not final_dir.exists():
        raise RuntimeError(f"MitoFinder final results directory not found: {final_dir}")

    gb_files = (
        list(final_dir.glob("*.gb")) +
        list(final_dir.glob("*.gbk")) +
        list(final_dir.glob("*.gbf"))
    )

    if not gb_files:
        raise RuntimeError(f"No GenBank file found in MitoFinder final results: {final_dir}")

    annotated_gb = gb_files[0]
    copied_gb = outdir / annotated_gb.name
    shutil.copy2(annotated_gb, copied_gb)

    return copied_gb
