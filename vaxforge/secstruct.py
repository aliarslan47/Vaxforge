"""GERÇEK ikincil yapı tahmini — S4PRED (PSIPRED grubu, tek-dizi model).

Homoloji/MSA kullanmadan yalnız birincil diziden C/H/E (coil/helix/strand)
tahmin eder (PyTorch ensemble, çevrimdışı). ss2 (PSIPRED VFORMAT) çıktısının
SS kolonundan %helix / %strand / %coil hesaplanır.

Atıf: Moffat L, Jones DT. Increasing the accuracy of single sequence prediction
methods using a deep semi-supervised learning framework. Bioinformatics.
2021;37(21):3744-3751. GitHub: psipred/s4pred.
"""

from __future__ import annotations

import functools
import subprocess
import sys
import tempfile
from pathlib import Path

_DIR = Path(__file__).resolve().parent.parent / "tools" / "s4pred"
_SCRIPT = _DIR / "run_model.py"
_WEIGHTS = _DIR / "weights" / "weights_1.pt"


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return _SCRIPT.exists() and _WEIGHTS.exists()


def _device() -> str:
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        return "cpu"


def predict(seq: str) -> dict:
    """Protein -> {helix_pct, strand_pct, coil_pct}. Araç yoksa boş sözlük."""
    if not available() or not seq:
        return {}
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(f">MEV\n{seq}\n")
        fpath = fh.name
    try:
        r = subprocess.run(
            [sys.executable, str(_SCRIPT), fpath, "-d", _device(), "-t", "ss2"],
            cwd=str(_DIR), capture_output=True, timeout=600, text=True,
        )
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return {}
    Path(fpath).unlink(missing_ok=True)
    # ss2: "pos AA SS p_coil p_helix p_strand" — SS kolonu (3.) C/H/E
    h = e = c = tot = 0
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit():
            ss = parts[2].upper()
            tot += 1
            if ss == "H":
                h += 1
            elif ss == "E":
                e += 1
            else:
                c += 1
    if not tot:
        return {}
    return {
        "helix_pct": round(100.0 * h / tot, 1),
        "strand_pct": round(100.0 * e / tot, 1),
        "coil_pct": round(100.0 * c / tot, 1),
    }
