"""GERÇEK sinyal peptidi tahmini — SignalP-5.0 (DTU, yerel).

tools/signalp-5.0b içindeki Go binary'yi çalıştırır. ÖNEMLİ: binary yalnızca
PATH'te 'signalp' basename'i ile çağrılınca doğru asset yolunu buluyor (go-bindata
quirk'i) + LD_LIBRARY_PATH lib/ (libtensorflow) gösterilmeli. Yoksa
available()=False -> funnel.py sinyal-peptidi heuristiğine düşer.
"""

from __future__ import annotations

import functools
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

_PKG = Path(__file__).resolve().parent.parent / "tools" / "signalp-5.0b"


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return (_PKG / "bin" / "signalp").exists() and (_PKG / "lib").exists()


def _org(profile: str) -> str:
    return {"bacteria": "gram-", "virus": "euk", "parasite": "euk"}.get(profile, "gram-")


def predict(proteins: list[tuple[str, str]], profile: str = "bacteria") -> dict[str, dict]:
    """[(id, seq)] -> {id: {sp_prob, has_sp}}. Çalışmazsa boş döner."""
    if not available() or not proteins:
        return {}
    idmap = {f"s{i}": pid for i, (pid, _) in enumerate(proteins)}
    with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
        for i, (_, seq) in enumerate(proteins):
            fh.write(f">s{i}\n{seq}\n")
        fpath = fh.name
    prefix = tempfile.mktemp(prefix="signalp_")
    env = dict(os.environ)
    env["PATH"] = str(_PKG / "bin") + os.pathsep + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = str(_PKG / "lib") + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    try:
        subprocess.run(["signalp", "-fasta", fpath, "-org", _org(profile),
                        "-format", "short", "-prefix", prefix],
                       capture_output=True, timeout=600, text=True, env=env)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return {}
    Path(fpath).unlink(missing_ok=True)
    out: dict[str, dict] = {}
    summary = Path(prefix + "_summary.signalp5")
    if summary.exists():
        for line in summary.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            if len(p) < 6:
                continue
            pid = idmap.get(p[0], p[0])
            try:
                sp_prob = float(p[2])
            except ValueError:
                sp_prob = 0.0
            out[pid] = {"sp_prob": round(sp_prob, 3), "has_sp": p[1] != "OTHER"}
        summary.unlink(missing_ok=True)
    return out
