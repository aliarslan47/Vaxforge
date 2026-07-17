"""GERÇEK toksisite tahmini — ToxinPred2 (Raghava lab, standalone).

tools/toxinpred2/toxinpred2.py yerel script'i ile çalışır. VARSAYILAN: Model 2
(Hybrid = AAC+RF ML + BLAST toksin-homolojisi + MERCI motif taraması) — yazarların
önerdiği, en yüksek doğruluklu model. Hybrid bağımlılıkları (blastp, perl, envfile,
Database, progs) yoksa Model 1'e (yalnız AAC+RF) düşer.

Model dosyası (RF_model) devasa olduğu için repoda değildir; ayrıca indirilir ve
bir kez kurulu sklearn'e migrate edilir (ifngamma.py sürüm-köprüsüyle aynı mantık).
Ayrıca upstream'de iki pandas-2.x uyumu (DataFrame.append, karışık-tip sum) düzeltildi.

Atıf: Sharma N, Naorem LD, Jain S, Raghava GPS. ToxinPred2: an improved method
for predicting toxicity of proteins. Brief Bioinform. 2022;23(5):bbac174.
"""

from __future__ import annotations

import csv
import functools
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_DIR = Path(__file__).resolve().parent.parent / "tools" / "toxinpred2"
_SCRIPT = _DIR / "toxinpred2.py"
_MODEL = _DIR / "RF_model"

# Model 2 (Hybrid) intermediate/çıktı dosyaları — çalıştırma sonrası temizlenir.
_JUNK = ["seq.aac", "seq.pred", "seq.out", "merci.txt", "merci_output.csv",
         "merci_hybrid.csv", "blast_hybrid.csv", "RES_1_6_6.out", "final_output",
         "Sequence_1"]


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return _SCRIPT.exists() and _MODEL.exists()


@functools.lru_cache(maxsize=1)
def hybrid_available() -> bool:
    """Model 2 (Hybrid) için BLAST + MERCI + veri dosyaları hazır mı."""
    env = _DIR / "envfile"
    return (shutil.which("blastp") is not None and shutil.which("perl") is not None
            and env.exists() and (_DIR / "Database" / "data.phr").exists()
            and (_DIR / "progs" / "MERCI_motif_locator.pl").exists()
            and (_DIR / "Database" / "pos_motif.txt").exists())


def _cleanup():
    for f in _JUNK:
        (_DIR / f).unlink(missing_ok=True)


def predict(peptides: list[str], threshold: float = 0.6) -> dict[str, dict]:
    """Peptit -> {ml_score, toxic, model}. Model 2 (Hybrid) tercih, yoksa Model 1.

    Not: script ara dosyalarını kendi dizinine yazdığı ve RF_model'i göreli
    yüklediği için cwd=_DIR ile çalıştırılır; ID'ler s{i} üzerinden eşlenir.
    """
    if not available() or not peptides:
        return {}
    peps = sorted({p for p in peptides if p})
    model = "2" if hybrid_available() else "1"
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fin:
        for i, p in enumerate(peps):
            fin.write(f">s{i}\n{p}\n")
        fpath = fin.name
    fout = fpath + ".csv"
    try:
        subprocess.run(
            [sys.executable, str(_SCRIPT), "-i", fpath, "-o", fout,
             "-m", model, "-d", "2", "-t", str(threshold)],
            cwd=str(_DIR), capture_output=True, timeout=1200, text=True,
        )
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        _cleanup()
        return {}
    out: dict[str, dict] = {}
    try:
        with open(fout, newline="") as fh:
            for row in csv.DictReader(fh):
                # Model 1: ID, Sequence, ML_Score, Prediction
                # Model 2: Subject, ML Score, MERCI Score, BLAST Score, Hybrid Score, Prediction
                rid = row.get("Subject") or row.get("ID")   # s{i}
                if not rid or not rid.startswith("s"):
                    continue
                try:
                    idx = int(rid[1:])
                except ValueError:
                    continue
                if idx >= len(peps):
                    continue
                score_s = (row.get("Hybrid Score") or row.get("ML_Score")
                           or row.get("ML Score") or "0")
                try:
                    score = float(score_s)
                except ValueError:
                    score = 0.0
                pred = (row.get("Prediction") or "").strip().lower()
                toxic = pred == "toxin" if pred else score >= threshold
                out[peps[idx]] = {"ml_score": round(score, 3), "toxic": toxic,
                                  "model": "hybrid" if model == "2" else "rf"}
    except FileNotFoundError:
        pass
    finally:
        Path(fpath).unlink(missing_ok=True)
        Path(fout).unlink(missing_ok=True)
        _cleanup()
    return out
