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

def run_cmd(cmd: List[str], log_file: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    if log_file:
        with open(log_file, "w", encoding="utf-8") as log:
            log.write("COMMAND:\n" + " ".join(cmd) + "\n\n")
            p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    else:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
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


def get_genetic_code(config: Dict[str, Any], default: int = 5) -> int:
    candidates = [
        safe_get(config, ["project", "genetic_code"], None),
        safe_get(config, ["genetic_code"], None),
        safe_get(config, ["annotation", "genetic_code"], None),
        safe_get(config, ["mitofinder", "organism_code"], None),
        safe_get(config, ["tools", "mitofinder", "organism_code"], None),
    ]
    for v in candidates:
        if v is not None and str(v).strip() != "":
            try:
                return int(v)
            except Exception:
                continue
    return int(default)
