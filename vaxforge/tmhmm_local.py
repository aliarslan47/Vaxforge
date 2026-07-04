"""Yerel TMHMM-2.0 (DTU) — transmembran heliks tahmini (opsiyonel, çevrimdışı).

Kullanıcı TMHMM-2.0c'yi (DTU) tools/ altına koyarsa gerçek TMHMM ile TM heliks
sayısı tahmin edilir. Yoksa available()=False -> funnel.py Kyte-Doolittle
hidropati yöntemine düşer (o da klasik/gerçek bir yöntem, TMHMM kadar hassas değil).

Not: DeepTMHMM yalnız bulutta (BioLib, hesap gerekir); tmhmm.py/pytmhmm ise
lisanslı model / derleme sorunları çıkardı — bu yüzden TMHMM-2.0 indirmesi
en temiz çevrimdışı gerçek yol.
"""

from __future__ import annotations

import functools
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent.parent / "tools"


def _wrapper() -> str | None:
    hits = sorted(_TOOLS.glob("tmhmm-*/bin/tmhmm")) + sorted(_TOOLS.glob("TMHMM-*/bin/tmhmm"))
    return str(hits[0]) if hits else None


@functools.lru_cache(maxsize=1)
def available() -> bool:
    w = _wrapper()
    return bool(w) and shutil.which("perl") is not None


def predict(proteins: list[tuple[str, str]]) -> dict[str, int]:
    """[(id, seq)] -> {id: tm_helis_sayısı}. Çalışmazsa boş döner."""
    if not available() or not proteins:
        return {}
    w = _wrapper()
    with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
        for pid, seq in proteins:
            fh.write(f">{pid}\n{seq}\n")
        fpath = fh.name
    workdir = tempfile.mkdtemp(prefix="tmhmm_")
    try:
        r = subprocess.run([w, "-short", fpath], capture_output=True, text=True,
                           timeout=600, cwd=workdir)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        shutil.rmtree(workdir, ignore_errors=True)
        return {}
    Path(fpath).unlink(missing_ok=True)
    shutil.rmtree(workdir, ignore_errors=True)
    out: dict[str, int] = {}
    # -short çıktısı: "<id>\tlen=..\tExpAA=..\tFirst60=..\tPredHel=K\tTopology=.."
    for line in r.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        pid = line.split()[0]
        m = re.search(r"PredHel=(\d+)", line)
        if m:
            out[pid] = int(m.group(1))
    return out
