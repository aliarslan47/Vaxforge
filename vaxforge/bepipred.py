"""GERÇEK lineer B-hücre epitop tahmini — BepiPred-1.0 (DTU, yerel).

BepiPred-1.0 (Larsen 2006): HMM + Parker hidrofilisite birleşimi. tcsh sarmalayıcı
+ statik binary; protein başına per-residue epitop skoru (GFF) üretir.

Kuruluysa (tools/bepipred-1.0/bepipred) gerçek skorları döndürür; yoksa
available()=False -> epitope.py klasik Kolaskar-Tongaonkar+Parker yöntemine düşer.
Varsayılan epitop eşiği 0.35 (BepiPred-1.0 standardı).
"""

from __future__ import annotations

import functools
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

_WRAPPER = Path(__file__).resolve().parent.parent / "tools" / "bepipred-1.0" / "bepipred"
DEFAULT_THRESHOLD = 0.35


@functools.lru_cache(maxsize=1)
def _configured() -> bool:
    """Sarmalayıcıdaki `setenv BP <yol>` GERÇEK bir dizini mi gösteriyor?

    NEDEN GEREKLİ: DTU paketleri kurulumda elle düzenlenmek üzere satıcı YER
    TUTUCU yoluyla gelir (`setenv BP /usr/cbs/packages/bepipred/1.0/...`).
    Düzenlenmezse sarmalayıcı dosyası VARDIR, tcsh de VARDIR, yani eski
    `available()` True döndürüyordu — ama araç `cannot find BP in ...` diyip
    çıkıyor, `predict_residues` BOŞ LİSTE dönüyor ve çağıran bunu "B-hücre
    epitopu bulunamadı" sanıyordu. Yani araç yokluğu, negatif BULGU gibi
    görünüyordu: projenin dürüstlük ilkesinin tam tersi.

    Aynı kusur netMHCpan/netMHCIIpan sarmalayıcılarında da vardı
    (bkz. netmhc_local._wrapper_home_ok).
    """
    if not _WRAPPER.exists():
        return False
    try:
        txt = _WRAPPER.read_text(errors="ignore")
    except Exception:
        return False
    m = re.search(r"^\s*setenv\s+BP\s+(\S+)", txt, re.M)
    if not m:
        return True          # desen yoksa karar verme (ileride paket değişirse)
    return Path(m.group(1)).is_dir()


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return _configured() and shutil.which("tcsh") is not None


def predict_residues(protein: str) -> list[float]:
    """Protein -> per-residue BepiPred epitop skorları. Hata olursa boş liste."""
    if not available() or len(protein) < 3:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
        fh.write(f">q\n{protein}\n")
        fpath = fh.name
    try:
        r = subprocess.run([str(_WRAPPER), fpath], capture_output=True,
                           text=True, timeout=300)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return []
    Path(fpath).unlink(missing_ok=True)
    scores: list[float] = []
    for line in r.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        # GFF: seqname source feature start end score ...
        if len(parts) >= 6:
            try:
                scores.append(float(parts[5]))
            except ValueError:
                continue
    return scores
