"""GERÇEK çözünürlük tahmini — Protein-Sol (Manchester, Warwicker grubu).

Birincil diziden 'scaled solubility' (0-1) tahmin eder; Niwa ve ark. 2009
hücre-serbest E. coli çözünürlük verisine fit. Eşik 0.45: üstü = ortalama E.
coli proteininden daha çözünür. Saf Perl pipeline, çevrimdışı.

Atıf: Hebditch M, Carballo-Amador MA, Charonis S, Curtis R, Warwicker J.
Protein-Sol: a web tool for predicting protein solubility from sequence.
Bioinformatics. 2017;33(19):3098-3100.
"""

from __future__ import annotations

import functools
import shutil
import subprocess
import uuid
from pathlib import Path

_DIR = (Path(__file__).resolve().parent.parent / "tools" / "protein-sol"
        / "protein-sol-sequence-prediction-software")
_WRAPPER = _DIR / "multiple_prediction_wrapper_export.sh"


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return _WRAPPER.exists() and (_DIR / "seq_reference_data.txt").exists() \
        and shutil.which("perl") is not None


def predict(seq: str) -> dict:
    """Protein -> {scaled_solubility, soluble, population_sol, percent_sol}.

    scaled_solubility: Protein-Sol ölçekli skoru (>0.45 → çözünür). Araç yoksa boş.
    Not: wrapper ara dosyaları cwd'ye (_DIR) yazar; benzersiz adla çalışıp temizleriz.
    """
    if not available() or not seq:
        return {}
    tag = f"vf_{uuid.uuid4().hex[:8]}"
    fin = _DIR / f"{tag}.fasta"
    fin.write_text(f">{tag}\n{seq}\n")
    generated = [fin, _DIR / f"{tag}.fasta_ORIGINAL", _DIR / "seq_prediction.txt",
                 _DIR / "seq_composition.txt", _DIR / "run.log",
                 _DIR / "seq_props.out", _DIR / "STYprops.out"]
    try:
        subprocess.run(["bash", str(_WRAPPER), fin.name],
                       cwd=str(_DIR), capture_output=True, timeout=300, text=True)
        pred = _DIR / "seq_prediction.txt"
        if not pred.exists():
            return {}
        scaled = pop = perc = None
        for line in pred.read_text().splitlines():
            if line.startswith("SEQUENCE PREDICTIONS,"):
                # SEQUENCE PREDICTIONS, ID, percent-sol, scaled-sol, population-sol, pI
                f = [x.strip() for x in line.split(",")]
                if len(f) >= 5:
                    perc, scaled, pop = float(f[2]), float(f[3]), float(f[4])
                break
        if scaled is None:
            return {}
        return {
            "scaled_solubility": round(scaled, 3),
            "soluble": scaled > 0.45,
            "population_sol": round(pop, 3) if pop is not None else None,
            "percent_sol": round(perc, 2) if perc is not None else None,
        }
    except Exception:
        return {}
    finally:
        for p in generated:
            p.unlink(missing_ok=True)
