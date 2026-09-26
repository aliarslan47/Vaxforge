"""ESMFold 3D yapı öngörüsü — ESMAtlas fold API (anahtarsız, GPU'suz).

Yerel ESMFold GPU + ~11GB model ister; bu makinede yok. Bunun yerine ESM
Metagenomic Atlas'ın herkese açık fold ucu kullanılır:
  POST https://api.esmatlas.com/foldSequence/v1/pdb/  (gövde = ham dizi, ≤400 aa)
Dönen PDB'de B-factor sütunu = pLDDT (kalıntı-başı güven, 0-100).

DÜRÜSTLÜK: ağ yoksa / dizi >400 aa / API hatası → None döner (uydurma yok);
çağıran 'ARAÇ_YOK' olarak dürüstçe raporlar. Bu bir ÖNGÖRÜDÜR — deneysel yapı değil.
"""

from __future__ import annotations

import socket
import urllib.request
from io import StringIO

API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_LEN = 400          # ESMAtlas fold ucu sınırı
_AA = set("ACDEFGHIKLMNPQRSTVWY")


def available(timeout: float = 3.0) -> bool:
    """api.esmatlas.com 443'e hızlı erişim var mı? (indirme/folding tetiklemez)."""
    try:
        with socket.create_connection(("api.esmatlas.com", 443), timeout=timeout):
            return True
    except OSError:
        return False


def _sanitize(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in _AA)


def fold(seq: str, timeout: int = 300) -> str | None:
    """Dizi → PDB metni (ESMAtlas). Hata/uzun/ağ-yok → None (dürüst)."""
    s = _sanitize(seq)
    if not s or len(s) > MAX_LEN:
        return None
    try:
        req = urllib.request.Request(API_URL, data=s.encode(), method="POST")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            pdb = resp.read().decode("utf-8", "replace")
        return pdb if "ATOM" in pdb else None
    except Exception:
        return None


def _ca_plddt(pdb: str) -> list[float]:
    """CA atomlarının B-factor'ı (pLDDT) — kalıntı-başı güven, 0-100 ölçeğinde.

    ESMAtlas fold ucu pLDDT'yi 0-1 ölçeğinde döndürebilir; 0-100'e normalize edilir
    (ESMFold standart ölçeği: >70 güvenilir, <50 düzensiz).
    """
    out: list[float] = []
    for line in pdb.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                out.append(float(line[60:66]))
            except ValueError:
                continue
    if out and max(out) <= 1.5:      # 0-1 ölçeği → 0-100
        out = [x * 100.0 for x in out]
    return out


def predict(seq: str, timeout: int = 300) -> dict | None:
    """Dizi → {pdb, n_res, plddt_mean, plddt_min, plddt_per_res}. Hata → None.

    pLDDT>70 güvenilir, >90 çok yüksek; <50 düzensiz/kararsız (ESMFold ölçeği).
    """
    pdb = fold(seq, timeout=timeout)
    if not pdb:
        return None
    plddt = _ca_plddt(pdb)
    if not plddt:
        return None
    return {
        "pdb": pdb,
        "n_res": len(plddt),
        "plddt_mean": round(sum(plddt) / len(plddt), 1),
        "plddt_min": round(min(plddt), 1),
        "plddt_per_res": [round(x, 1) for x in plddt],
    }
