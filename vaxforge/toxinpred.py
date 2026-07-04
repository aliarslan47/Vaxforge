"""GERÇEK toksisite tahmini — ToxinPred2 (Raghava lab, pip paketi).

`toxinpred2` CLI kuruluysa gerçek RF (AAC) modeliyle tahmin eder; yoksa
available()=False döner ve survival.py proxy'ye düşer.
Model 1 (AAC+RF): BLAST/MERCI gerektirmez, hızlı, çevrimdışı.
"""

from __future__ import annotations

import csv
import functools
import shutil
import subprocess
import tempfile
from pathlib import Path


@functools.lru_cache(maxsize=1)
def _cli() -> str | None:
    return shutil.which("toxinpred2")


def available() -> bool:
    return _cli() is not None


def predict(peptides: list[str], threshold: float = 0.6) -> dict[str, dict]:
    """Peptit -> {ml_score, toxic}. ToxinPred2 model 1. Hata/araç yoksa boş döner."""
    cli = _cli()
    if not cli or not peptides:
        return {}
    peps = sorted({p for p in peptides if p})
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fin:
        for i, p in enumerate(peps):
            fin.write(f">s{i}\n{p}\n")
        fpath = fin.name
    fout = fpath + ".csv"
    try:
        subprocess.run([cli, "-i", fpath, "-o", fout, "-m", "1", "-d", "2"],
                       capture_output=True, timeout=600, text=True)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return {}
    out: dict[str, dict] = {}
    try:
        with open(fout, newline="") as fh:
            for row in csv.DictReader(fh):
                seq = row.get("Sequence") or row.get("seq")
                if not seq:
                    continue
                try:
                    score = float(row.get("ML_Score", "0"))
                except ValueError:
                    score = 0.0
                pred = (row.get("Prediction") or "").strip().lower()
                toxic = pred == "toxin" if pred else score >= threshold
                out[seq] = {"ml_score": round(score, 3), "toxic": toxic}
    except FileNotFoundError:
        pass
    finally:
        Path(fpath).unlink(missing_ok=True)
        Path(fout).unlink(missing_ok=True)
    return out
