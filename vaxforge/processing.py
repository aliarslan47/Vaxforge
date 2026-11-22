"""Antijen işleme omurgası seçimi (proteozomal kesim + TAP).

NEDEN BU KATMAN: işleme omurgası değiştirilebilir ve her omurganın skor ÖLÇEĞİ
kendine özgüdür (ör. IEDB `proteasome_score` log-ölçeklidir ve değerlerin
%68'i 1'i aşar; 0-1 olasılık varsayan bir normalizasyondan geçirmek kesim
bileşeninin çözünürlüğünün üçte ikisini kırpar — ölçüm: iedb_processing modül
başlığı). Bu yüzden normalizasyonu omurganın KENDİSİ yapar: çağıran yer modülü
alır ve `processing_norm`/`METHOD`'u ondan okur.

Şu an tek omurga IEDB processing API'dir. Yoksa None → MHC-I işleme metrikleri
dürüstçe ölçülmemiş kalır (ölçülmemiş ≠ olumsuz bulgu).
"""

from __future__ import annotations

from types import ModuleType

from . import iedb_processing


def backend() -> ModuleType | None:
    """Bu koşuda kullanılacak omurga modülü ya da None.

    Dönen modül şu üçünü sunar: `METHOD`, `predict(proteins)`,
    `processing_norm(cle, tap)`.
    """
    if iedb_processing.available():
        return iedb_processing
    return None


def available() -> bool:
    return backend() is not None


def method() -> str:
    b = backend()
    return b.METHOD if b else "ARAÇ_YOK (kesim+TAP ölçülmedi)"
