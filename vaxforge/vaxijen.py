"""VaxiJen antijenite çapraz-kontrolü — bağımsız 2. antijenite kanıtı.

Ana antijenite IApred'dir (Miles 2025; VaxiJen'den daha iyi raporlanır). Bu modül
İKİNCİ, BAĞIMSIZ bir kanıt sağlar: VaxiJen'in yayınlanmış ÖZELLİK yöntemi (Sandberg
z-ölçekleri + auto-cross-covariance + lojistik) — antigen_acc üzerinden yerel olarak.
İki bağımsız sinyalin uyuşması adayın antijenitesine güveni artırır.

Resmi VaxiJen sunucusu (ddg-pharmfac.net) bot-korumalı (HTTP 403) ve otomasyona
kapalı; opsiyonel canlı sorgu 403-toleranslıdır → başarısızlıkta None (dürüst ARAÇ_YOK).
DÜRÜSTLÜK: yerel skor VaxiJen'in YÖNTEMİDİR, resmi eğitilmiş SVM'i DEĞİL.
"""

from __future__ import annotations

import urllib.parse
import urllib.request

from . import antigen_acc

THRESHOLD = 0.5   # antigen_acc lojistik olasılık eşiği


def available() -> bool:
    """Yerel VaxiJen-yöntemi (antigen_acc) modeli mevcut mu?"""
    return antigen_acc.available()


def local_score(seq: str) -> dict:
    """VaxiJen-tarzı ACC skoru (yerel). {score, antigenic, method}."""
    s = antigen_acc.predict(seq)
    return {
        "score": s,
        "antigenic": s >= THRESHOLD,
        "method": "GERÇEK (VaxiJen-tarzı ACC — 2. kanıt; resmi SVM değil)",
    }


def server_score(seq: str, target: str = "bacteria", timeout: int = 25) -> float | None:
    """Resmi VaxiJen sunucusuna opsiyonel sorgu. 403/hata → None (dürüst)."""
    try:
        data = urllib.parse.urlencode({"seq": seq, "target": target}).encode()
        req = urllib.request.Request(
            "http://www.ddg-pharmfac.net/vaxijen/scripts/VaxiJen_scripts/VaxiJen3.pl",
            data=data, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            html = resp.read().decode("utf-8", "replace")
        import re
        m = re.search(r"Overall Prediction.*?=\s*([0-9.]+)", html)
        return float(m.group(1)) if m else None
    except Exception:
        return None


def cross_check(seq: str, try_server: bool = False) -> dict:
    """Bağımsız 2. antijenite kanıtı. {score, antigenic, method, server(opsiyonel)}."""
    out = local_score(seq)
    if try_server:
        srv = server_score(seq)
        out["server_score"] = srv
        out["server_note"] = ("resmi VaxiJen erişilemedi (403/hata)" if srv is None
                               else "resmi VaxiJen sunucusu")
    return out
