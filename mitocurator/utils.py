from __future__ import annotations
from pathlib import Path
import subprocess
import shutil
import yaml
from typing import Any, Dict, List, Optional

def load_config(path: str | Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)

def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

def run_cmd(cmd, log_file=None, check=True, cwd=None):
    import subprocess
    from pathlib import Path

    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)

        with open(log_file, "w") as log:
            p = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=cwd,
            )
    else:
        p = subprocess.run(cmd, cwd=cwd)

    if check and p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(cmd)}")

    return p

def which_or_path(value: str | None) -> str | None:
    if not value:
        return None
    if "/" in value:
        return value if Path(value).exists() else None
    return shutil.which(value)

def safe_get(d: Dict[str, Any], keys: List[str], default=None):
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur
