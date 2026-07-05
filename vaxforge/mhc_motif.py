"""MHC bağlanma yarığı cep/anchor motifi — YALNIZ YORUMLAMA (sıralamayı değiştirmez).

Bir peptidin sunulabilmesi için anchor kalıntılarının MHC yarığının ceplerine
(MHC-I: P2 → B cebi, C-terminal PΩ → F cebi) oturması gerekir. Bu bilgi NetMHCpan'in
%rank'ında zaten örtük olarak vardır (allel, yarığı astarlayan 34 kalıntının
pseudo-sekansıyla temsil edilir); bu modül onu YALNIZ AÇIKLAMA amacıyla görünür kılar.

allele_anchor_motif(): allelin anchor pozisyonlarında hangi amino asitleri tercih
ettiğini EZBERDEN DEĞİL, yerel NetMHCpan'e tek-substitüsyon taraması yaptırarak
AMPİRİK çıkarır (veri-güdümlü, tekrarlanabilir, çevrimdışı, NetMHCpan'e atıflanabilir
— [[reverse-vaccinology-pipeline]] 'LLM literatür madenciliği YOK' kararıyla uyumlu).
Sonuç diske önbelleklenir.

Atıflar: anchor kavramı Falk & Rammensee 1991; cep adlandırması Saper-Bjorkman 1991;
küratörlü motif DB SYFPEITHI (Rammensee 1999); pseudo-sekans Nielsen 2007.
"""

from __future__ import annotations

import json
from pathlib import Path

from . import netmhc_local

AA = "ACDEFGHIKLMNPQRSTVWY"
_BG = "AAALAAAAV"          # hafif bağlanan 9-mer arka plan (tarama için)
_CACHE = Path(__file__).resolve().parent.parent / "outputs" / ".mhc_motif_cache.json"

METHOD = "AMPİRİK (NetMHCpan tek-substitüsyon taraması)"


def anchor_positions(mhc_class: str, peplen: int) -> list[int]:
    """1-tabanlı primer anchor pozisyonları. MHC-I: P2 + C-terminal PΩ."""
    if mhc_class == "mhc_i":
        return [2, peplen] if peplen >= 3 else []
    # MHC-II çekirdek anchorları (P1/P4/P6/P9) çekirdek hizalaması gerektirir;
    # NetMHCIIpan çekirdeği metrik olarak tutulmuyor -> ampirik motif verilmez.
    return []


def peptide_anchors(pep: str, mhc_class: str) -> dict[str, str]:
    """Peptidin kendi anchor kalıntıları, ör. {'P2':'L', 'PΩ':'V'}."""
    pep = (pep or "").upper()
    pos = anchor_positions(mhc_class, len(pep))
    if not pos:
        return {}
    out: dict[str, str] = {}
    for p in pos:
        label = "PΩ" if p == len(pep) else f"P{p}"
        out[label] = pep[p - 1]
    return out


def _load_cache() -> dict:
    try:
        return json.loads(_CACHE.read_text())
    except Exception:
        return {}


def _save_cache(c: dict) -> None:
    try:
        _CACHE.parent.mkdir(parents=True, exist_ok=True)
        _CACHE.write_text(json.dumps(c, ensure_ascii=False))
    except Exception:
        pass


def allele_anchor_motif(allele: str, mhc_class: str = "mhc_i",
                        top: int = 4) -> dict[str, list[str]] | None:
    """Allelin anchor pozisyonlarında tercih ettiği amino asitler (ampirik).

    Yerel NetMHCpan'e arka-plan 9-mer üzerinde P2 ve PΩ'yi tek tek 20 amino aside
    değiştirip tarar; en iyi bağlanmayı veren kalıntıları döndürür. Yalnız MHC-I;
    araç yoksa None. Sonuç diske önbelleklenir.
    """
    if mhc_class != "mhc_i" or not allele:
        return None
    cache = _load_cache()
    key = f"{allele}|{mhc_class}"
    if key in cache:
        return cache[key]
    if not netmhc_local.runnable("mhc_i"):
        return None

    L = len(_BG)
    scan_pos = [2, L]                       # P2, PΩ
    peps, meta = [], []
    for pos in scan_pos:
        for aa in AA:
            p = list(_BG)
            p[pos - 1] = aa
            peps.append("".join(p))
            meta.append((pos, aa))
    res = netmhc_local.predict(peps, [allele], "mhc_i", 100.0)
    if not res:
        return None
    by_pos: dict[int, list[tuple[str, float]]] = {p: [] for p in scan_pos}
    for (pos, aa), pep in zip(meta, peps):
        d = res.get(pep)
        if d:
            by_pos[pos].append((aa, d["rank"]))
    motif: dict[str, list[str]] = {}
    for pos in scan_pos:
        ranked = sorted(by_pos[pos], key=lambda x: x[1])
        label = "PΩ" if pos == L else f"P{pos}"
        motif[label] = [aa for aa, _ in ranked[:top]]
    cache[key] = motif
    _save_cache(cache)
    return motif


def annotate(peptide, allele: str, mhc_class: str) -> None:
    """Peptide.metrics'e anchor kalıntıları + allel cep tercihi + uyum notunu yazar.

    SALT YORUMLAMA: candidacy/sıralamayı etkilemez.
    """
    anchors = peptide_anchors(peptide.seq, mhc_class)
    if not anchors:
        return
    peptide.metrics["anchor_residues"] = anchors
    motif = allele_anchor_motif(allele, mhc_class)
    if motif:
        peptide.metrics["allele_anchor_motif"] = motif
        # her anchor için peptidin kalıntısı allelin tercih listesinde mi?
        matches = [lbl for lbl, aa in anchors.items() if aa in motif.get(lbl, [])]
        peptide.metrics["anchor_match"] = f"{len(matches)}/{len(anchors)}"
        pretty = "; ".join(f"{lbl}={aa}" for lbl, aa in anchors.items())
        pref = "; ".join(f"{lbl}∈{{{','.join(motif[lbl])}}}" for lbl in motif)
        peptide.metrics["anchor_note"] = (
            f"Peptit anchorları: {pretty}. {allele} cep tercihi: {pref}. "
            f"Uyum: {len(matches)}/{len(anchors)} (NetMHCpan %rank cep-uyumunu zaten puanlar)."
        )
        peptide.methods["anchor_motif"] = METHOD
    else:
        # ampirik motif yok (araç yok / MHC-II) -> yalnız peptit anchorları + kavram
        pretty = "; ".join(f"{lbl}={aa}" for lbl, aa in anchors.items())
        peptide.metrics["anchor_note"] = (
            f"Peptit anchorları: {pretty}. Anchor kalıntıları MHC yarığının "
            f"ceplerine oturmalıdır (Falk-Rammensee 1991); bağlanma uyumu NetMHCpan "
            f"%rank'ında değerlendirilir."
        )
        peptide.methods["anchor_motif"] = "kavramsal (peptit anchorları; allel motifi hesaplanmadı)"
