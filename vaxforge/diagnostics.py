"""Tanılama kaydı — yakalanan istisnaların İZ BIRAKMASINI sağlar.

NEDEN: pipeline aşamaları `try/except` ile sarılıdır (ROADMAP kural 2) ki tek bir
aracın çökmesi koşuyu düşürmesin. Ancak istisna yalnız `type(e).__name__` olarak
not edilirse iki durum birbirinden AYIRT EDİLEMEZ:

  (a) araç kurulu değil  → dürüst ARAÇ_YOK, beklenen davranış
  (b) araç kurulu ama çöktü → gerçek hata, sessizce kayboluyor

`record()` her yakalanan istisnayı `meta["diagnostics"]` listesine yapısal olarak
yazar (aşama, istisna tipi, mesaj, son traceback karesi) ve aynı anda `logging` ile
bildirir. `report.py` tüm meta'yı run.json'a serileştirdiği için kayıt otomatik
olarak koşu çıktısına ve provenance zincirine girer.

DÜRÜSTLÜK: burada hiçbir şey uydurulmaz; yalnızca gerçekten oluşan istisna yazılır.
"""

from __future__ import annotations

import logging
import os
import traceback

LOGGER_NAME = "vaxforge"

log = logging.getLogger(LOGGER_NAME)


def setup(level: str | int | None = None) -> logging.Logger:
    """Kök 'vaxforge' logger'ını bir kez yapılandırır (idempotent).

    Seviye: `level` argümanı, yoksa VAXFORGE_LOG ortam değişkeni, yoksa INFO.
    Streamlit/FastAPI kendi handler'larını kurduğu için propagate kapatılmaz;
    yalnız handler yoksa bir StreamHandler eklenir (çift satır olmaz).
    """
    lvl = level or os.environ.get("VAXFORGE_LOG") or "INFO"
    if isinstance(lvl, str):
        lvl = getattr(logging, lvl.upper(), logging.INFO)
    log.setLevel(lvl)
    if not log.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter(
            "%(asctime)s %(levelname)s [%(name)s] %(message)s",
            datefmt="%H:%M:%S"))
        log.addHandler(h)
    return log


def _last_frame(exc: BaseException) -> str:
    """İstisnanın en son traceback karesi: 'dosya:satır in fonksiyon'."""
    tb = traceback.extract_tb(exc.__traceback__)
    if not tb:
        return ""
    f = tb[-1]
    return f"{os.path.basename(f.filename)}:{f.lineno} in {f.name}"


def record(meta: dict | None, stage: str, exc: BaseException,
           expected: bool = False) -> dict:
    """Yakalanan istisnayı meta['diagnostics']'e yazar + loglar; kaydı döndürür.

    `expected=True` yalnızca aracın yokluğu gibi BEKLENEN düşüşler için kullanılır
    (seviye: INFO). Beklenmeyen istisna WARNING seviyesinde bildirilir, çünkü
    koşu devam etse bile bir araç gerçekten çökmüş demektir.
    """
    entry = {
        "stage": stage,
        "error": type(exc).__name__,
        "message": str(exc)[:500],
        "where": _last_frame(exc),
        "expected": bool(expected),
    }
    if isinstance(meta, dict):
        meta.setdefault("diagnostics", []).append(entry)
    msg = (f"{stage}: {entry['error']}: {entry['message']}"
           f"{' @ ' + entry['where'] if entry['where'] else ''}")
    log.info(msg) if expected else log.warning(msg)
    return entry


def summarize(meta: dict | None) -> dict:
    """Rapor/UI için kısa özet: kaç beklenmeyen hata, hangi aşamalarda."""
    items = (meta or {}).get("diagnostics") or []
    unexpected = [d for d in items if not d.get("expected")]
    return {
        "total": len(items),
        "unexpected": len(unexpected),
        "stages": sorted({str(d.get("stage")) for d in unexpected}),
    }
