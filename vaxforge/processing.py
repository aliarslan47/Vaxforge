"""Antijen işleme omurgası seçimi: NetCTL-1.2 (yerel) ya da IEDB processing API.

NEDEN BU KATMAN: iki tahminci aynı iki bileşeni (proteozomal kesim + TAP)
üretir ama ÖLÇEKLERİ farklıdır (NetCTL cle = 0-1 olasılık, IEDB
proteasome_score = log-ölçekli, %68'i 1'i aşıyor). Normalizasyonu hangi omurga
koştuysa onun yapması gerekir: IEDB skorlarını netctl.processing_norm'dan
geçirmek kesim bileşeninin çözünürlüğünün üçte ikisini kırpıyor (ölçüm:
iedb_processing modül başlığı). Bu yüzden çağıran yer omurganın KENDİSİNİ alır
ve `processing_norm`/`METHOD`'u ondan okur.

SIRA: NetCTL yerelde GERÇEKTEN koşuyorsa (duman testi geçerse) o tercih edilir —
makalenin özgün yöntemi, çevrimdışı ve hızlı. Koşmuyorsa IEDB API'ye düşülür.
İkisi de yoksa None → MHC-I işleme metrikleri dürüstçe ölçülmemiş kalır.
"""

from __future__ import annotations

from types import ModuleType

from . import iedb_processing, netctl


def backend() -> ModuleType | None:
    """Bu koşuda kullanılacak omurga modülü ya da None.

    Dönen modül şu üçünü sunar: `METHOD`, `predict(proteins)`,
    `processing_norm(cle, tap)`.
    """
    if netctl.available():
        return netctl
    if iedb_processing.available():
        return iedb_processing
    return None


def available() -> bool:
    return backend() is not None


def method() -> str:
    b = backend()
    return b.METHOD if b else "ARAÇ_YOK (kesim+TAP ölçülmedi)"
