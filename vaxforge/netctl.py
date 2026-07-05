"""Adım 4 (MHC-I sonrası) — antijen işleme: proteozomal C-terminal kesim + TAP.

GERÇEK araç: NetCTL-1.2 (Larsen ve ark. 2005, Eur J Immunol 35:2295). Bir
peptidin sunulabilmesi için yalnız MHC-I'e bağlanması yetmez; önce proteozom
tarafından doğru C-terminalden kesilmeli ve TAP ile ER'ye taşınmalıdır. NetCTL
bunun için iki allel-BAĞIMSIZ skor üretir:

  * cle — C-terminal proteozomal kesim olasılığı (0-1, sinir ağı).
  * tap — TAP taşıma yakınlığı (log-odds, ~ -3..+3).

TASARIM KARARI: NetCTL'in KENDİ eski MHC bağlanma tahmincisini (yalnız 12
supertype) KULLANMIYORUZ — MHC bağlanmayı NetMHCpan-4.2 yapıyor. Buradan yalnız
allel-bağımsız `cle` + `tap` alınır ve funnel-survivor proteinlerinin her 9-mer'i
için C-terminal pozisyona göre indekslenir; epitope._mhc bunları MHC-I
peptitlerimize (peptidin C-terminaline denk gelen 9-mer) eşler.

NetCTL yalnız 9-mer tarar. 8/10/11-mer peptitlerimiz için, aynı C-terminale sahip
9-mer'in cle/tap'ı kullanılır: kesim genuinely bir C-terminal sinyalidir (peptit
uzunluğundan bağımsız geçerli), TAP için ise yaklaşıklıktır (dürüstçe böyle
etiketlenir).
"""

from __future__ import annotations

import functools
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent.parent / "tools"

METHOD = "GERÇEK (NetCTL-1.2 kesim+TAP, Larsen 2005)"

# NetCTL'in içsel göreli ağırlıkları: kesim (wc=0.15) : TAP (wt=0.05) = 3:1.
# MHC bileşenini (we=1.0) DIŞARIDA bırakırız (onu NetMHCpan yapıyor).
_W_CLE, _W_TAP = 0.15, 0.05

_ROW = re.compile(
    r"\bID\s+(?P<id>\S+)\s+pep\s+(?P<pep>\S+).*?\bcle\s+(?P<cle>-?\d+\.\d+)"
    r"\s+tap\s+(?P<tap>-?\d+\.\d+)"
)


def _wrapper() -> str | None:
    hits = sorted(_TOOLS.glob("netCTL-*/netCTL")) + sorted(_TOOLS.glob("netctl-*/netCTL"))
    return str(hits[0]) if hits else None


@functools.lru_cache(maxsize=1)
def available() -> bool:
    w = _wrapper()
    return bool(w) and shutil.which("tcsh") is not None


def processing_norm(cle: float | None, tap: float | None) -> float:
    """cle (0-1) + tap (log-odds) -> 0-1 antijen işleme skoru.

    TAP log-odds sigmoid ile 0-1'e çekilir; NetCTL'in kesim:TAP=3:1 göreli
    ağırlığıyla birleştirilir. İkisi de yoksa 0.0.
    """
    if cle is None and tap is None:
        return 0.0
    parts, wsum = 0.0, 0.0
    if cle is not None:
        parts += _W_CLE * max(0.0, min(1.0, cle))
        wsum += _W_CLE
    if tap is not None:
        parts += _W_TAP * (1.0 / (1.0 + math.exp(-tap)))
        wsum += _W_TAP
    return round(parts / wsum, 4) if wsum else 0.0


def predict(proteins: list[tuple[str, str]]) -> dict[str, dict[int, dict[str, float]]]:
    """[(id, seq)] -> {protein_id: {cterm_pos: {"cle":.., "tap":..}}}.

    cterm_pos = 9-mer'in 0-tabanlı C-terminal kalıntı indeksi (start+8). Araç
    yoksa / hata olursa boş döner (çağıran yer sessizce atlar).
    """
    if not available() or not proteins:
        return {}
    w = _wrapper()
    # Sentetik kısa ID (NetCTL uzun/özel ID'leri kısaltabilir) -> gerçek ID eşlemesi.
    id_map = {f"s{i}": pid for i, (pid, _seq) in enumerate(proteins)}
    with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
        for i, (_pid, seq) in enumerate(proteins):
            fh.write(f">s{i}\n{seq}\n")
        fpath = fh.name
    try:
        r = subprocess.run([w, fpath], capture_output=True, text=True, timeout=900)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return {}
    finally:
        Path(fpath).unlink(missing_ok=True)

    out: dict[str, dict[int, dict[str, float]]] = {}
    # Her protein için satırlar sırayla gelir; pozisyon = o protein içindeki
    # 9-mer sırası (0-tabanlı start). C-terminal = start + 8.
    order: dict[str, int] = {}
    for line in r.stdout.splitlines():
        m = _ROW.search(line)
        if not m:
            continue
        sid = m.group("id")
        pid = id_map.get(sid, sid)
        start = order.get(sid, 0)
        order[sid] = start + 1
        cterm = start + 8
        out.setdefault(pid, {})[cterm] = {
            "cle": round(float(m.group("cle")), 4),
            "tap": round(float(m.group("tap")), 4),
        }
    return out
