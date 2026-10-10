"""Adım 6b — Adaylık puanı.

Her hayatta kalan peptidin metriklerini config'teki ağırlıklarla tek skora
indirger. Bir peptit tüm bileşenlere sahip olmayabilir (ör. B-hücre peptidinin
MHC skoru yoktur); o durumda mevcut bileşenler üzerinden ağırlıklar yeniden
normalize edilir. (Yapısal kalite/pLDDT bileşeni, yapısal doğrulama adımı
şimdilik çıkarıldığı için kullanılmıyor.)
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
    if m.get("virulence") is not None:  # kaynak proteinin virülans kanıtı (VFDB skoru)
        comp["virulence"] = float(m["virulence"])
    # Viral zarf/yüzey glikoproteini sinyali (yalnız virüs profilinde set edilir; SignalP+TM imzası).
    # Bakteride hiç set edilmez → comp'ta yer almaz → bakteri skoru DEĞİŞMEZ.
    if m.get("surface_glyco") is not None:
        comp["surface_exposure"] = float(m["surface_glyco"])
    if p.kind == "MHC-I":
        comp["mhc_i_binding"] = float(m.get("mhc_score", 0))
        if m.get("immunogenicity_norm") is not None:
            comp["immunogenicity"] = float(m["immunogenicity_norm"])
        if m.get("processing_norm") is not None:
            comp["antigen_processing"] = float(m["processing_norm"])
    if p.kind == "MHC-II":
        comp["mhc_ii_binding"] = float(m.get("mhc_score", 0))
        if m.get("ifn_gamma_norm") is not None:
            comp["ifn_gamma"] = float(m["ifn_gamma_norm"])
    if p.kind == "B":
        comp["antigenicity"] = max(comp.get("antigenicity", 0), float(m.get("bcell_score", 0)))
    if m.get("coverage_frac") is not None:
        comp["organism_coverage"] = float(m["coverage_frac"])
    # Suş konservasyonu (yalnız kullanıcı ek suş verdiyse hesaplanır → varsa bileşen).
    # Epitop-penceresi korunmuşluğu (0-1); korunmuş epitop = geniş koruma → yüksek istenir.
    if m.get("epitope_conservation") is not None:
        comp["conservation"] = float(m["epitope_conservation"])
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
        # Moleküler mimikri: TAM eşleşme artık pipeline'da eleniyor (sert kapı).
        # Buraya düşen tam-eşleşme yalnız kapı kapalıyken (discard_exact_match
        # false) gelir; o durumda eski yumuşak ceza korunur.
        if p.metrics.get("self_mimicry"):
            val *= 0.9
            p.notes.append("moleküler mimikri (konak self tam eşleşme) → skor ×0.9")
        # BLOSUM-yakın mimikri: elenmez, konfigürasyondan okunan ceza uygulanır
        # (TCR dejenerasyonu — birebir olmayan peptit de çapraz-reaktif olabilir).
        mim_factor = p.metrics.get("mimicry_penalty_factor")
        if mim_factor:
            val *= float(mim_factor)
            p.notes.append(f"kısmi mimikri (BLOSUM62 yakın) → skor ×{mim_factor}")
        # Protein düzeyi konak-homolojisi: kaynak antijen konağa benziyorsa
        # (protein-bütünü kimlik ≥ min_identity) ceza. Eleme YOK — ölçüldü:
        # protein-bütünü kimlik epitop-düzeyi riskin vekili değil (enolaz %49.8
        # kimlik, 9-mer'lerinin %1.9'u self). Risk mimikri kapısında kapanıyor.
        hh_factor = p.metrics.get("host_homology_penalty_factor")
        if hh_factor:
            val *= float(hh_factor)
            p.notes.append(f"kaynak antijen konak-benzeri işaretli → skor ×{hh_factor}")
        # Toksisite: aynı mantıkta yumuşak ceza. Katsayı survival.run'da
        # config'ten (toxicity.toxicity_penalty) okunup metriğe yazılır; epitop
        # düzeyinde ELEME YOK, bkz. survival.py modül açıklaması.
        tox_factor = p.metrics.get("toxicity_penalty_factor")
        if tox_factor:
            val *= float(tox_factor)
            p.notes.append(f"toksik işaretli (ToxinPred3) → skor ×{tox_factor}")
        p.candidacy = round(val, 4)
        p.metrics["score_components"] = comp
    peptides.sort(key=lambda x: x.candidacy, reverse=True)
    return peptides
