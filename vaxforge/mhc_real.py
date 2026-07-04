"""Gerçek MHC-I bağlanma tahmini — MHCflurry (eğitilmiş sinir ağı, GPU'suz).

MHCflurry kuruluysa temsili insan HLA-I panelinde GERÇEK presentation tahmini
yapar. Kurulu değilse available()=False döner ve epitope.py proxy'ye düşer.

Not: MHCflurry insan HLA-I içindir. Diğer türler (fare H-2, sığır BoLA…) için
NetMHCpan gerekir; onlar bu modülde 'çalıştırılmadı' sayılır. Böylece çok-tür
kapsamı için ŞU AN insan-içi HLA kapsamı gerçek, tür-arası kısım eksik/etiketli.
"""

from __future__ import annotations

import functools

# IEDB referans temsili insan HLA-I paneli (supertype kapsayıcı)
HLA_PANEL = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*40:01",
]


@functools.lru_cache(maxsize=1)
def _predictor():
    try:
        from mhcflurry import Class1PresentationPredictor
        return Class1PresentationPredictor.load()
    except Exception:
        return None


def available() -> bool:
    return _predictor() is not None


def supported_alleles() -> set:
    pred = _predictor()
    if pred is None:
        return set()
    try:
        return set(pred.supported_alleles)
    except Exception:
        return set()


def predict(peptides: list[str], alleles: list[str] | None = None,
            weak_percentile: float = 2.0) -> dict[str, dict]:
    """Peptit -> {best_percentile, best_allele, n_alleles} verilen HLA panelinde.

    alleles verilmezse varsayılan insan HLA_PANEL kullanılır. MHCflurry yalnız
    insan HLA-I desteklediğinden, desteklenmeyen alleller sessizce atlanır.
    n_alleles: kaç allelin bu peptidi sunduğu (percentile ≤ eşik).
    """
    pred = _predictor()
    if pred is None or not peptides:
        return {}
    panel = list(alleles) if alleles else HLA_PANEL
    supported = supported_alleles()
    panel = [a for a in panel if a in supported]
    if not panel:
        return {}
    peps = sorted({p for p in peptides if 8 <= len(p) <= 15})
    if not peps:
        return {}
    df = pred.predict(peptides=peps, alleles={a: [a] for a in panel}, verbose=0)
    out: dict[str, dict] = {}
    for pep, grp in df.groupby("peptide"):
        best_row = grp.loc[grp["presentation_percentile"].idxmin()]
        n_present = int((grp["presentation_percentile"] <= weak_percentile).sum())
        out[pep] = {
            "best_percentile": round(float(best_row["presentation_percentile"]), 3),
            "best_allele": str(best_row["sample_name"]),
            "presentation_score": round(float(best_row["presentation_score"]), 4),
            "n_alleles": n_present,
            "panel_size": len(panel),
        }
    return out
