"""Adım 6b — Adaylık puanı.

Her hayatta kalan peptidin metriklerini config'teki ağırlıklarla tek skora
indirger. Bir peptit tüm bileşenlere sahip olmayabilir (ör. B-hücre peptidinin
MHC skoru yoktur); o durumda mevcut bileşenler üzerinden ağırlıklar yeniden
normalize edilir. Yapısal kalite (pLDDT) ağır/deferred adımdan gelir; yoksa
bileşen dışında bırakılır.
"""

from __future__ import annotations

from .models import Peptide


def _components(p: Peptide) -> dict[str, float]:
    """Peptidin normalize (0-1) skor bileşenleri."""
    m = p.metrics
    comp: dict[str, float] = {}
    parent_ag = p.metrics.get("parent_antigenicity")
    if parent_ag is not None:
        comp["antigenicity"] = float(parent_ag)
    if p.kind == "MHC-I":
        comp["mhc_i_binding"] = float(m.get("mhc_score", 0))
    if p.kind == "MHC-II":
        comp["mhc_ii_binding"] = float(m.get("mhc_score", 0))
    if p.kind == "B":
        comp["antigenicity"] = max(comp.get("antigenicity", 0), float(m.get("bcell_score", 0)))
    if m.get("conservation") is not None:
        comp["conservation"] = float(m["conservation"])
    if m.get("coverage_frac") is not None:
        comp["organism_coverage"] = float(m["coverage_frac"])
    if m.get("plddt") is not None:
        comp["structural_quality"] = float(m["plddt"]) / 100
    return comp


def _dedupe(peptides: list[Peptide]) -> list[Peptide]:
    """Aynı (dizi, tip) peptidini tekilleştir (ilk görülen kalır)."""
    seen, out = set(), []
    for p in peptides:
        key = (p.seq, p.kind)
        if key not in seen:
            seen.add(key)
            out.append(p)
    return out


def score(peptides: list[Peptide], weights: dict[str, float]) -> list[Peptide]:
    peptides = _dedupe(peptides)
    for p in peptides:
        comp = _components(p)
        if not comp:
            p.candidacy = 0.0
            continue
        wsum = sum(weights.get(k, 0) for k in comp) or 1.0
        val = sum(weights.get(k, 0) * v for k, v in comp.items()) / wsum
        p.candidacy = round(val, 4)
        p.metrics["score_components"] = comp
    peptides.sort(key=lambda x: x.candidacy, reverse=True)
    return peptides
