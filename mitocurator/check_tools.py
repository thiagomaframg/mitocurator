from __future__ import annotations
from pathlib import Path
from .utils import which_or_path, run_cmd, safe_get, ensure_dir

BASIC_TOOLS = ["minimap2", "samtools", "bcftools", "blastn", "tblastn", "blastp", "mafft", "seqkit", "arwen"]

def check_tools(config: dict, outdir: Path):
    ensure_dir(outdir)
    rows = []
    tools = config.get("tools", {}) or {}

    for t in BASIC_TOOLS:
        configured = tools.get(t, t)
        path = which_or_path(configured)
        status = "OK" if path else "MISSING"
        rows.append([t, status, configured, path or "."])

    mf = tools.get("mitofinder", {}) or {}
    if mf.get("enabled", False):
        mode = mf.get("mode", "python_interpreter")
        status = "UNKNOWN"
        detail = "."
        command = "."
        try:
            if mode == "python_interpreter":
                py2 = mf.get("python2")
                script = mf.get("script")
                py2_ok = Path(py2).exists() if py2 else False
                script_ok = Path(script).exists() if script else False
                command = f"{py2} {script} --help"
                if py2_ok and script_ok:
                    p = run_cmd([py2, script, "--help"], check=False)
                    status = "OK" if p.returncode in (0, 1, 2) else "CHECK"
                    detail = f"python2={py2}; script={script}; returncode={p.returncode}"
                else:
                    status = "MISSING"
                    detail = f"python2_exists={py2_ok}; script_exists={script_ok}"
            elif mode == "conda_env":
                conda = mf.get("conda_executable", "conda")
                env = mf.get("conda_env")
                script = mf.get("script")
                conda_path = which_or_path(conda)
                command = f"{conda} run -n {env} {script} --help"
                if conda_path and env and script and Path(script).exists():
                    p = run_cmd([conda_path, "run", "-n", env, script, "--help"], check=False)
                    status = "OK" if p.returncode in (0, 1, 2) else "CHECK"
                    detail = f"conda={conda_path}; env={env}; script={script}; returncode={p.returncode}"
                else:
                    status = "MISSING"
                    detail = f"conda={conda_path}; env={env}; script={script}"
            elif mode == "wrapper":
                wrapper = mf.get("wrapper")
                command = f"{wrapper} --help"
                if wrapper and Path(wrapper).exists():
                    p = run_cmd([wrapper, "--help"], check=False)
                    status = "OK" if p.returncode in (0, 1, 2) else "CHECK"
                    detail = f"wrapper={wrapper}; returncode={p.returncode}"
                else:
                    status = "MISSING"
                    detail = f"wrapper={wrapper}"
            else:
                status = "ERROR"
                detail = f"Unsupported mode: {mode}"
        except Exception as e:
            status = "ERROR"
            detail = str(e)
        rows.append(["mitofinder", status, mode, detail])

    illum_polish = (safe_get(config, ["polish", "illumina"], None) or {})
    if illum_polish:
        bwa2 = which_or_path("bwa-mem2")
        bwa  = which_or_path("bwa")
        if bwa2:
            rows.append(["bwa-mem2", "OK", "bwa-mem2", bwa2])
        elif bwa:
            rows.append(["bwa-mem2", "MISSING", "bwa-mem2", "."])
            rows.append(["bwa", "OK (fallback)", "bwa", bwa])
        else:
            rows.append(["bwa-mem2", "MISSING", "bwa-mem2", "."])
            rows.append(["bwa",      "MISSING",  "bwa",      "."])

    outfile = outdir / "tool_check.tsv"
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("tool\tstatus\tconfigured\tresolved_or_detail\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    return outfile
