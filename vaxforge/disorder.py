"""GERÇEK içsel düzensizlik (intrinsic disorder) tahmini — metapredict V3.

Holehouse lab'ın konsensüs düzensizlik tahmincisi (PyTorch, çevrimdışı). Kayıt
gerektirmez (IUPred3 aksine). Per-rezidü skor + aracın kendi IDR (düzensiz bölge)
çağrılarını verir; %disordered bu bölgelerin kapsamından hesaplanır.

Atıf: Emenecker RJ, Griffith D, Holehouse AS. Metapredict: a fast, accurate, and
easy-to-use predictor of consensus disorder and structure. Biophys J.
2021;120(20):4312-4319. (V2/V3: Emenecker ve ark., bioRxiv 2022).
"""

from __future__ import annotations

import functools


@functools.lru_cache(maxsize=1)
def available() -> bool:
    try:
        import metapredict  # noqa: F401
        return True
    except Exception:
        return False


def predict(seq: str) -> dict:
    """Protein -> {percent_disordered, mean_score, n_idrs, idr_boundaries}.

    %disordered: metapredict'in kendi düzensiz-bölge (IDR) çağrılarının toplam
    uzunluğu / dizi uzunluğu. Araç yoksa boş sözlük döner.
    """
    if not available() or not seq:
        return {}
    import warnings

    import numpy as np
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import metapredict as meta
        r = meta.predict_disorder(seq, return_domains=True)
    n = len(seq)
    idrs = list(getattr(r, "disordered_domain_boundaries", []) or [])
    covered = sum(max(0, int(e) - int(s)) for s, e in idrs)
    scores = np.asarray(getattr(r, "disorder", []), dtype=float)
    return {
        "percent_disordered": round(100.0 * covered / n, 1) if n else 0.0,
        "mean_score": round(float(scores.mean()), 3) if scores.size else 0.0,
        "n_idrs": len(idrs),
        "idr_boundaries": [[int(s), int(e)] for s, e in idrs],
    }
