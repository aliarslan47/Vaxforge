"""MHC bağlanma yarığı cep/anchor motifi — YALNIZ YORUMLAMA (sıralamayı değiştirmez).

Bir peptidin sunulabilmesi için anchor kalıntılarının MHC yarığının ceplerine
oturması gerekir. Bu bilgi NetMHCpan'in %rank'ında zaten örtük olarak vardır
(allel, yarığı astarlayan 34 kalıntının pseudo-sekansıyla temsil edilir); bu
modül onu YALNIZ AÇIKLAMA amacıyla görünür kılar.

ANCHOR POZİSYONLARI ALLELE ÖZGÜDÜR (2026-11-22 düzeltmesi)
----------------------------------------------------------
Bu modülün ilk sürümü anchor pozisyonlarını TÜM alleller için P2 + C-terminal
PΩ olarak SABİTLİYORDU. Bu biyolojik olarak yanlıştır: hangi pozisyonun anchor
olduğunu, yarığı astarlayan polimorfik kalıntıların oluşturduğu cepler belirler
(P2 -> B cebi, PΩ -> F cebi, ama ayrıca P3 -> D cebi ve P5 -> E cebi). Cepler
allelden allele değişir; konak yalnız HANGİ ALLEL REPERTUARINA sahip olduğu
için önemlidir. Aynı konakta bile alleller farklı pozisyonları kullanır.

Bu 16 allelde yerel NetMHCpan ile ÖLÇÜLDÜ (pozisyon başına tek-substitüsyon
taraması, üç nötr arka plan, log10 %rank dinamik aralığı):

    HLA-A*02:01   P2(1.42) P9(1.28)              -> eski varsayım doğru
    HLA-A*01:01   P9(2.27) P3(1.13) P2(0.88)     -> P3 > P2, eski kod P3'ü kaçırıyordu
    HLA-B*08:01   P9(1.04) P5(1.00) P3(0.78)     -> P2 ilk üçte bile değil
    H-2-Kb        P9(1.48) P3(1.10) P5(0.87)     -> P2 anchor DEĞİL
    H-2-Db        P9(1.63) P5(1.46)              -> klasik Db P5 anchor'ı
    SLA-1*0101    P9(1.80) P3(1.65) P2(1.18)     -> P3 ~ P9
    BoLA-3:00101  P9(2.29) P2(1.17) P3(1.13)     -> P3 ~ P2

16 allelin 8'inde P2+PΩ dışında güçlü bir anchor pozisyonu bulundu; 6'sında P2
ilk iki belirleyici pozisyonda değildi. Dört konağın hepsi etkilendi.

Bu yüzden anchor pozisyonları artık EZBERDEN DEĞİL, her allel için AMPİRİK
olarak ölçülür. Araç yoksa P2+PΩ'ye düşülür ama bu durum raporda "varsayılan
(ölçülmedi)" olarak DÜRÜSTÇE etiketlenir.

Atıflar: anchor kavramı Falk & Rammensee 1991; cep adlandırması Saper-Bjorkman
1991; küratörlü motif DB SYFPEITHI (Rammensee 1999); pseudo-sekans Nielsen 2007.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from . import netmhc_local

AA = "ACDEFGHIKLMNPQRSTVWY"

# Tarama arka planları. Tek bir arka plan bir anchor'ı maskeleyebildiği için
# ÜÇ nötr dizi kullanılır ve dinamik aralıklar ortalanır.
_BG_I = ("AAALAAAAV", "SGAGLGAGS", "KTAGQGTYA")            # 9-mer
_BG_II = ("AAAAAKAAAAVAAAA", "SGAGLGAGSAGLSGA", "KTAGQGTYAQGTKAG")  # 15-mer

_CACHE_VERSION = 2   # v1 = sabit P2+PΩ; v2 = ölçülen pozisyonlar
_CACHE = Path(__file__).resolve().parent.parent / "outputs" / ".mhc_motif_cache.json"

# Bir pozisyonun anchor sayılması için dinamik aralık ölçütleri.
# spread = log10(en kötü %rank) - log10(en iyi %rank); anchor'da büyük olur.
ANCHOR_SPREAD_MIN = 0.50    # mutlak alt sınır (yaklaşık 3 kat rank oynaması)
ANCHOR_SPREAD_FRAC = 0.45   # en güçlü pozisyonun bu kesrinden büyük olmalı
MAX_ANCHORS = 4             # en fazla bu kadar pozisyon anchor olarak bildirilir

METHOD_MEASURED = "AMPİRİK (NetMHCpan pozisyon-başına tek-substitüsyon taraması)"
METHOD_DEFAULT = "VARSAYILAN (P2+PΩ; allel-özgü ölçüm YAPILMADI — araç yok)"
METHOD = METHOD_MEASURED  # geriye uyumluluk


def _lg(rank) -> float:
    return math.log10(max(float(rank), 0.01))


def _backgrounds(mhc_class: str) -> tuple[str, ...]:
    return _BG_I if mhc_class == "mhc_i" else _BG_II


# --------------------------------------------------------------------------- #
# Önbellek
# --------------------------------------------------------------------------- #
def _load_cache() -> dict:
    try:
        c = json.loads(_CACHE.read_text())
    except Exception:
        return {"_version": _CACHE_VERSION}
    if c.get("_version") != _CACHE_VERSION:
        # Sürüm atladı -> eski (sabit pozisyonlu) girdiler GEÇERSİZ, atılır.
        return {"_version": _CACHE_VERSION}
    return c


def _save_cache(c: dict) -> None:
    try:
        c["_version"] = _CACHE_VERSION
        _CACHE.parent.mkdir(parents=True, exist_ok=True)
        _CACHE.write_text(json.dumps(c, ensure_ascii=False))
    except Exception:
        pass


# --------------------------------------------------------------------------- #
# Ampirik pozisyon taraması
# --------------------------------------------------------------------------- #
def scan_positions(allele: str, mhc_class: str = "mhc_i") -> dict | None:
    """Her pozisyonun bağlanmaya katkısını ölçer (allel-özgü anchor keşfi).

    Arka plan peptidinde pozisyon p'deki kalıntıyı 20 amino aside tek tek
    değiştirip %rank ölçer. Anchor pozisyonunda rank büyük ölçüde oynar; anchor
    olmayanda neredeyse sabit kalır.

    Dönen: {"L": <tarama uzunluğu>, "pos": {p: {"spread": x, "top": [aa...]}}}
    Araç yoksa None.
    """
    if not allele or not netmhc_local.runnable(mhc_class):
        return None
    bgs = _backgrounds(mhc_class)
    L = len(bgs[0])
    peps: list[str] = []
    meta: list[tuple[str, int, str]] = []
    for bg in bgs:
        for pos in range(1, L + 1):
            for aa in AA:
                p = list(bg)
                p[pos - 1] = aa
                peps.append("".join(p))
                meta.append((bg, pos, aa))
    res = netmhc_local.predict(peps, [allele], mhc_class, 100.0)
    if not res:
        return None

    acc: dict[int, dict[str, list[tuple[str, float]]]] = {}
    for (bg, pos, aa), pep in zip(meta, peps, strict=True):
        d = res.get(pep)
        if d:
            acc.setdefault(pos, {}).setdefault(bg, []).append((aa, d["rank"]))

    out: dict[int, dict] = {}
    for pos, bg_map in acc.items():
        spreads: list[float] = []
        prefs: dict[str, list[float]] = {}
        for lst in bg_map.values():
            if len(lst) < 10:          # yetersiz veri -> bu arka planı atla
                continue
            ranks = [r for _, r in lst]
            spreads.append(_lg(max(ranks)) - _lg(min(ranks)))
            for aa, r in lst:
                prefs.setdefault(aa, []).append(_lg(r))
        if not spreads:
            continue
        ranked = sorted(((a, sum(v) / len(v)) for a, v in prefs.items()),
                        key=lambda x: x[1])
        out[pos] = {
            "spread": round(sum(spreads) / len(spreads), 3),
            "top": [a for a, _ in ranked[:4]],
        }
    if not out:
        return None
    return {"L": L, "pos": out}


def _label(pos: int, scan_len: int) -> str:
    """Tarama pozisyonunu etikete çevirir. Son pozisyon C-terminal = PΩ."""
    return "PΩ" if pos == scan_len else f"P{pos}"


def allele_profile(allele: str, mhc_class: str = "mhc_i") -> dict | None:
    """Allelin ölçülmüş anchor profili (önbellekli).

    Dönen: {"measured": True, "scan_len": L, "anchors": {label: {...}},
            "all": {label: spread}}  ya da araç yoksa None.
    """
    if not allele:
        return None
    cache = _load_cache()
    key = f"{allele}|{mhc_class}"
    if key in cache:
        return cache[key]
    scan = scan_positions(allele, mhc_class)
    if scan is None:
        return None
    L, pos = scan["L"], scan["pos"]
    best = max(v["spread"] for v in pos.values())
    cutoff = max(ANCHOR_SPREAD_MIN, best * ANCHOR_SPREAD_FRAC)
    chosen = sorted((p for p, v in pos.items() if v["spread"] >= cutoff),
                    key=lambda p: -pos[p]["spread"])[:MAX_ANCHORS]
    prof = {
        "measured": True,
        "scan_len": L,
        "anchors": {
            _label(p, L): {"pos": p, "spread": pos[p]["spread"],
                           "prefers": pos[p]["top"]}
            for p in sorted(chosen)
        },
        "all": {_label(p, L): v["spread"] for p, v in sorted(pos.items())},
    }
    cache[key] = prof
    _save_cache(cache)
    return prof


# --------------------------------------------------------------------------- #
# Pozisyon/kalıntı erişimi
# --------------------------------------------------------------------------- #
def anchor_positions(mhc_class: str, peplen: int,
                     allele: str | None = None) -> list[int]:
    """Bu allel ve peptit uzunluğu için 1-tabanlı anchor pozisyonları.

    Allel verilir ve ölçüm mümkünse ÖLÇÜLEN pozisyonlar; aksi hâlde MHC-I için
    varsayılan P2 + PΩ (dürüstçe 'varsayılan' etiketlenir), MHC-II için boş.
    """
    if peplen < 3:
        return []
    if allele:
        prof = allele_profile(allele, mhc_class)
        if prof:
            out: list[int] = []
            for lbl, info in prof["anchors"].items():
                p = peplen if lbl == "PΩ" else info["pos"]
                if 1 <= p <= peplen and p not in out:
                    out.append(p)
            return sorted(out)
    if mhc_class == "mhc_i":
        return sorted({2, peplen})
    # MHC-II: çekirdek hizalaması bilinmeden varsayılan pozisyon verilmez.
    return []


def peptide_anchors(pep: str, mhc_class: str,
                    allele: str | None = None) -> dict[str, str]:
    """Peptidin anchor kalıntıları, ör. {'P2':'L', 'P5':'F', 'PΩ':'V'}."""
    pep = (pep or "").upper()
    pos = anchor_positions(mhc_class, len(pep), allele)
    if not pos:
        return {}
    return {("PΩ" if p == len(pep) else f"P{p}"): pep[p - 1] for p in pos}


def allele_anchor_motif(allele: str, mhc_class: str = "mhc_i",
                        top: int = 4) -> dict[str, list[str]] | None:
    """Allelin anchor pozisyonlarında tercih ettiği amino asitler (ampirik).

    Geriye uyumlu imza: {label: [aa, ...]} döner. Araç yoksa None.
    """
    prof = allele_profile(allele, mhc_class)
    if not prof:
        return None
    return {lbl: info["prefers"][:top] for lbl, info in prof["anchors"].items()}


# --------------------------------------------------------------------------- #
# Peptide'e yazma
# --------------------------------------------------------------------------- #
def annotate(peptide, allele: str, mhc_class: str, host: str | None = None) -> None:
    """Peptide.metrics'e anchor kalıntıları + allel cep tercihi + dayanağı yazar.

    SALT YORUMLAMA: candidacy/sıralamayı etkilemez.

    `host` verilirse, anchor'ların HANGİ TÜRÜN hangi alleline göre belirlendiği
    `anchor_basis` alanına yazılır. Anchor'lar allele özgü olduğu için bu bilgi
    olmadan tablo yorumlanamaz; rapor bu alanı tablonun altına basar.
    """
    prof = allele_profile(allele, mhc_class)
    anchors = peptide_anchors(peptide.seq, mhc_class, allele)
    if not anchors:
        return
    peptide.metrics["anchor_residues"] = anchors

    host_txt = f"{host} / " if host else ""
    if prof:
        motif = {lbl: info["prefers"] for lbl, info in prof["anchors"].items()}
        peptide.metrics["allele_anchor_motif"] = motif
        matches = [lbl for lbl, aa in anchors.items() if aa in motif.get(lbl, [])]
        peptide.metrics["anchor_match"] = f"{len(matches)}/{len(anchors)}"
        spreads = "; ".join(
            f"{lbl}={info['spread']}" for lbl, info in prof["anchors"].items())
        peptide.metrics["anchor_positions"] = ",".join(prof["anchors"].keys())
        peptide.metrics["anchor_basis"] = (
            f"{host_txt}{allele} ({mhc_class.upper().replace('MHC_','MHC-')}); "
            f"anchor pozisyonları bu allel için ÖLÇÜLDÜ "
            f"(dinamik aralık log10 %rank: {spreads})"
        )
        pretty = "; ".join(f"{lbl}={aa}" for lbl, aa in anchors.items())
        pref = "; ".join(f"{lbl}∈{{{','.join(motif[lbl])}}}" for lbl in motif)
        peptide.metrics["anchor_note"] = (
            f"Peptit anchorları: {pretty}. {allele} cep tercihi: {pref}. "
            f"Uyum: {len(matches)}/{len(anchors)} "
            f"(NetMHCpan %rank cep-uyumunu zaten puanlar)."
        )
        peptide.methods["anchor_motif"] = METHOD_MEASURED
    else:
        pretty = "; ".join(f"{lbl}={aa}" for lbl, aa in anchors.items())
        peptide.metrics["anchor_positions"] = ",".join(anchors.keys())
        peptide.metrics["anchor_basis"] = (
            f"{host_txt}{allele}; anchor pozisyonları ÖLÇÜLMEDİ — "
            f"varsayılan P2+PΩ kullanıldı. Anchor pozisyonları ALLELE ÖZGÜDÜR; "
            f"bu allelde gerçek anchor başka pozisyonlarda olabilir."
        )
        peptide.metrics["anchor_note"] = (
            f"Peptit anchorları (varsayılan P2+PΩ): {pretty}. Anchor kalıntıları "
            f"MHC yarığının ceplerine oturmalıdır (Falk-Rammensee 1991); "
            f"bağlanma uyumu NetMHCpan %rank'ında değerlendirilir."
        )
        peptide.methods["anchor_motif"] = METHOD_DEFAULT


def basis_lines(peptides) -> list[str]:
    """Tablonun altına basılacak 'hangi tür/allel' satırları (tekilleştirilmiş).

    Anchor'lar allele özgü olduğundan, bir anchor tablosu hangi konağın hangi
    alleline göre hesaplandığı yazılmadan yorumlanamaz.
    """
    seen: dict[str, None] = {}
    for p in peptides or []:
        b = (getattr(p, "metrics", {}) or {}).get("anchor_basis")
        if b:
            seen.setdefault(str(b), None)
    return list(seen)
