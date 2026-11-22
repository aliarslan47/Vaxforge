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


# ── pLDDT güven yayılımı ───────────────────────────────────────────────────
# NEDEN (2026-11-22): BoHV-1 MEV'inin pLDDT ortalaması 55.6 ölçüldü, yani
# "güvenilir" eşiğinin (70) ALTINDA. Rapor bu sayıyı diğer satırlar gibi
# yazıyordu; oysa YAPIDAN TÜRETİLEN her sonuç (RSA/yüzey yüzdesi, konformasyonel
# B-hücre yamaları) o yapının güvenine BAĞLIDIR. Düşük pLDDT'li bir modelden
# okunan "yüzey epitobu" ölçülmüş bir bulgu gibi görünüyordu.
#
# ÖNEMLİ BAĞLAM: MEV, linkerlarla birleştirilmiş yapay bir kimeradır; parçaları
# doğal bir katlanma birimi DEĞİL. Böyle bir dizide düşük pLDDT BEKLENEN
# durumdur, bir kusur değil. Bu yüzden rapor "kötü" demez: ilgili alt sonuçları
# "düşük güven" diye işaretler ve nedenini yazar.

#: ESMFold ölçeği: >70 güvenilir, 50-70 düşük, <50 düzensiz.
PLDDT_CONFIDENT = 70.0
PLDDT_DISORDERED = 50.0


def confidence(plddt_mean: float | None) -> dict:
    """Ortalama pLDDT → {level, reliable, label_tr, label_en}."""
    if plddt_mean is None:
        return {"level": "unknown", "reliable": False,
                "label_tr": "ölçülmedi", "label_en": "not measured"}
    v = float(plddt_mean)
    if v >= 90:
        return {"level": "very_high", "reliable": True,
                "label_tr": "çok yüksek güven", "label_en": "very high confidence"}
    if v >= PLDDT_CONFIDENT:
        return {"level": "confident", "reliable": True,
                "label_tr": "güvenilir", "label_en": "confident"}
    if v >= PLDDT_DISORDERED:
        return {"level": "low", "reliable": False,
                "label_tr": "DÜŞÜK güven (eşik 70'in altında)",
                "label_en": "LOW confidence (below the 70 threshold)"}
    return {"level": "disordered", "reliable": False,
            "label_tr": "ÇOK DÜŞÜK — düzensiz model",
            "label_en": "VERY LOW — disordered model"}


def region_plddt(plddt_per_res: list[float] | None,
                 start: int, end: int) -> float | None:
    """[start, end] (1-TABANLI, iki uç dahil) aralığının ortalama pLDDT'si."""
    if not plddt_per_res:
        return None
    a = max(0, int(start) - 1)
    b = min(len(plddt_per_res), int(end))
    if b <= a:
        return None
    win = plddt_per_res[a:b]
    return round(sum(win) / len(win), 1)
