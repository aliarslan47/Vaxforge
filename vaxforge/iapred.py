"""GERÇEK antijenite — IApred (Miles ve ark., ImmunoInformatics 2025).

Açık kaynak, host-ve-patojen bağımsız intrinsik antijenite tahmincisi; SVM,
918 yüksek-antijenik proteinle eğitilmiş (bakteri/virüs/mantar/protozoa/helmint).
VaxiJen 2.0/3.0 ve ANTIGENpro'yu geçiyor (ROC AUC 0.761). VaxiJen otomatikleşmediği
için (web 403) bunu kullanıyoruz. GitHub: sebamiles/IAPred.

tools/IAPred/IApred.py çalıştırılır; skor 'Intrinsic_Antigenicity_Score'
(≈ -3..+3; >0.3 yüksek, <-0.3 düşük). Yoksa available()=False -> antigen_acc/heuristik.
"""

from __future__ import annotations

import csv
import functools
import math
import subprocess
import sys
import tempfile
from pathlib import Path

_DIR = Path(__file__).resolve().parent.parent / "tools" / "IAPred"
_SCRIPT = _DIR / "IApred.py"


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return _SCRIPT.exists() and (_DIR / "models").exists()


def _norm(raw: float) -> float:
    """IApred ham skorunu (≈-3..3) 0-1 aralığına (sigmoid) çevirir."""
    return round(1 / (1 + math.exp(-raw)), 3)


def predict(proteins: list[tuple[str, str]]) -> dict[str, dict]:
    """[(id, seq)] -> {id: {score(raw), norm(0-1), category, antigenic}}."""
    if not available() or not proteins:
        return {}
    idmap = {f"s{i}": pid for i, (pid, _) in enumerate(proteins)}
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
        for i, (_, seq) in enumerate(proteins):
            fh.write(f">s{i}\n{seq}\n")
        fpath = fh.name
    out_csv = fpath + ".csv"
    try:
        subprocess.run([sys.executable, str(_SCRIPT), fpath, out_csv],
                       cwd=str(_DIR), capture_output=True, timeout=900, text=True)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return {}
    Path(fpath).unlink(missing_ok=True)
    out: dict[str, dict] = {}
    p = Path(out_csv)
    if p.exists():
        with p.open() as fh:
            for row in csv.DictReader(fh):
                pid = idmap.get(row.get("Header", ""), row.get("Header", ""))
                try:
                    raw = float(row.get("Intrinsic_Antigenicity_Score", "0"))
                except ValueError:
                    raw = 0.0
                out[pid] = {"score": raw, "norm": _norm(raw),
                            "category": row.get("Antigenicity_Category", ""),
                            "antigenic": raw > 0.3}
        p.unlink(missing_ok=True)
    return out
