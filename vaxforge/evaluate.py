"""Aday-başına TÜM araç sonuçları + eşiğe göre GEÇTİ/GEÇEMEDİ değerlendirmesi.

Her aday peptit için, uygulanabilir her aracın çıktısını (kaynak proteinin funnel
sonuçları + peptit-seviyesi tahminler) tek bir satır listesine döker. Eşik DEĞERLERİ
config'ten gelir (meta['thresholds']); karşılaştırma YÖNÜ ('≥ geçer' / '≤ geçer')
her metriğin sabit bilimsel özelliğidir ve burada kodlanır.

Durum anlamı:
  - Sert filtre (hard): geçemeyen aday zaten ELENİR; sağ kalanlar 'GEÇTİ' görünür
    (şeffaflık için gösterilir).
  - Skor aracı (soft): geçememek ELEMEZ, yalnız adaylık puanını etkiler → aday
    listede kalsa da satırda 'GEÇEMEDİ' görülebilir.
  - Yorum (—): eşiği olmayan bilgilendirici satır (ör. TAP, anchor motifi).
"""

from __future__ import annotations

PASS, FAIL, NA = "✅ GEÇTİ", "❌ GEÇEMEDİ", "—"


def thr_lookup(meta: dict) -> dict:
    """meta['thresholds'] -> {(tool,param): (value, unit, hard_filter)}."""
    out = {}
    for r in meta.get("thresholds", []):
        out[(r.get("tool"), r.get("param"))] = (
            r.get("value"), r.get("unit", ""), bool(r.get("hard_filter")))
    return out


def _row(tool, value, cutoff, status, method, hard=False, kind="skor"):
    return {"tool": tool, "value": value, "cutoff": cutoff,
            "status": status, "method": method, "hard": hard, "kind": kind}


def report_subset(peptides, top_n: int = 15) -> list[tuple]:
    """Yayın raporuna (PDF/HTML) girecek aday alt kümesi.

    Tüm adayları uzun uzadıya dökmek raporu şişirir; bunun yerine:
      - en iyi `top_n` aday (adaylık puanına göre sıralı), +
      - literatürde (IEDB) bilinen bir epitopla eşleşen HER aday (top_n dışında olsa da).
    Tam liste ayrıca candidates_full.xlsx'e yazılır.

    Dönüş: [(orijinal_sıra, peptide, sebep)] — sebep ∈ {'top', 'literatür'};
    orijinal sıra numarası korunur (tam listedeki yeri).
    """
    out, seen = [], set()
    for rank, p in enumerate(peptides, 1):
        if rank <= top_n:
            out.append((rank, p, "top"))
            seen.add(id(p))
    for rank, p in enumerate(peptides, 1):
        if id(p) in seen:
            continue
        if (p.metrics.get("iedb") or {}).get("matched") is True:
            out.append((rank, p, "literatür"))
            seen.add(id(p))
    return out


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def candidate_rows(p, thr: dict) -> list[dict]:
    """Bir aday için sıralı araç-sonuç satırları (pipeline sırasıyla)."""
    m = p.metrics
    rows: list[dict] = []

    def gv(tool, param, default=None):
        return thr.get((tool, param), (default, "", False))[0]

    # ---- Kaynak protein: funnel araçları -----------------------------------
    loc = m.get("localization")
    if loc is not None:
        allowed = gv("localization", "allowed") or []
        st = PASS if loc in allowed else (FAIL if allowed else NA)
        rows.append(_row("Hücre-altı lokalizasyon (DeepLoc)", loc,
                         f"∈ {{{', '.join(allowed)}}}" if allowed else "—",
                         st if allowed else NA, m.get("method_localization", ""), hard=True))
    if m.get("signalp") is not None:
        v = _num(m.get("signalp")); c = _num(gv("signal_peptide", "min_probability"))
        st = (PASS if v >= c else FAIL) if (v is not None and c is not None) else NA
        rows.append(_row("Sinyal peptidi (SignalP)", round(v, 3) if v is not None else v,
                         f"≥ {c}", st, m.get("method_signalp", "")))
    if m.get("tm_helices") is not None:
        v = _num(m.get("tm_helices")); c = _num(gv("transmembrane", "max_helices"))
        st = (PASS if v <= c else FAIL) if (v is not None and c is not None) else NA
        rows.append(_row("Transmembran heliks (TMHMM)", int(v) if v is not None else v,
                         f"≤ {int(c) if c is not None else c}", st,
                         m.get("method_tm", ""), hard=True))
    hh = m.get("human_homology")
    if hh is not None:
        if hh == "not_run":
            rows.append(_row("İnsan homolojisi (diamond)", "çalıştırılmadı (DB yok)",
                             f"< %{gv('human_homology', 'min_identity')}", NA,
                             "ZORUNLU — DB bağlı değil", hard=True))
        else:
            v = _num(hh); c = _num(gv("human_homology", "min_identity"))
            st = (PASS if v < c else FAIL) if (v is not None and c is not None) else NA
            rows.append(_row("İnsan homolojisi (diamond, otoimmünite)",
                             f"%{round(v, 1) if v is not None else v}", f"< %{c}", st,
                             m.get("method_human_homology", ""), hard=True))
    ag = _num(m.get("parent_antigenicity"))
    if ag is not None:
        c = _num(gv("antigenicity_vaxijen", "threshold"))
        st = (PASS if ag >= c else FAIL) if c is not None else NA
        cat = m.get("antigenicity_category")
        rows.append(_row("Antijenite — kaynak protein (IApred)",
                         f"{round(ag, 3)}" + (f" ({cat})" if cat else ""),
                         f"≥ {c}", st, m.get("method_antigenicity", "")))

    # ---- Peptit-seviyesi ---------------------------------------------------
    if p.kind == "MHC-I":
        rank = _num(m.get("pseudo_rank"))
        weak = _num(gv("mhc_class_i", "rank_weak")); strong = _num(gv("mhc_class_i", "rank_strong"))
        if rank is not None and weak is not None:
            st = PASS if rank <= weak else FAIL
            lvl = " (güçlü)" if (strong is not None and rank <= strong) else (" (zayıf)" if rank <= weak else "")
            rows.append(_row("MHC-I bağlanma %rank (NetMHCpan)", f"{rank}{lvl}",
                             f"≤ {weak} (güçlü ≤ {strong})", st, m.get("mhc_score") or p.methods.get("mhc_score", "")))
        if m.get("immunogenicity") is not None:
            v = _num(m.get("immunogenicity"))
            st = PASS if (v is not None and v > 0) else FAIL
            rows.append(_row("CD8+ immünojenite (IEDB Calis 2013)", v, "> 0", st,
                             p.methods.get("immunogenicity", "")))
        if m.get("cleavage") is not None:
            v = _num(m.get("cleavage")); c = _num(gv("netctl", "min_cleavage"))
            st = (PASS if v >= c else FAIL) if c is not None else NA
            rows.append(_row("Proteozomal kesim cle (NetCTL)", v, f"≥ {c}", st,
                             p.methods.get("processing", "")))
            rows.append(_row("TAP taşıma (NetCTL)", m.get("tap"), "bilgi", NA,
                             p.methods.get("processing", "")))
        if m.get("anchor_note"):
            rows.append(_row("MHC yarığı anchor/cep motifi", m.get("anchor_match", ""),
                             "yorum", NA, p.methods.get("anchor_motif", "")))
    elif p.kind == "MHC-II":
        rank = _num(m.get("pseudo_rank"))
        weak = _num(gv("mhc_class_ii", "rank_weak")); strong = _num(gv("mhc_class_ii", "rank_strong"))
        if rank is not None and weak is not None:
            st = PASS if rank <= weak else FAIL
            rows.append(_row("MHC-II bağlanma %rank (NetMHCIIpan)", rank,
                             f"≤ {weak} (güçlü ≤ {strong})", st, p.methods.get("mhc_score", "")))
    elif p.kind == "B":
        v = _num(m.get("bcell_score")); c = _num(gv("bcell_epitope", "min_score"))
        st = (PASS if v >= c else FAIL) if (v is not None and c is not None) else NA
        rows.append(_row("B-hücre epitop skoru (BepiPred)", v, f"≥ {c}", st,
                         p.methods.get("bcell_score", "")))

    # konak/tür kapsamı (bilgilendirici)
    if m.get("host_coverage") is not None:
        rows.append(_row("Konak/tür kapsamı (sunan konak sayısı)",
                         f"{m.get('host_coverage')} konak", "bilgi", NA, ""))

    # ---- Sağ kalım: sert filtreler -----------------------------------------
    if "allergen" in m:
        is_alg = bool(m.get("allergen"))
        rows.append(_row("Alerjenite (FAO/WHO 6-mer)", "ALERJEN" if is_alg else "temiz",
                         "non-alerjen", FAIL if is_alg else PASS,
                         p.methods.get("allergen", ""), hard=True))
    if m.get("toxicity") is not None:
        v = _num(m.get("toxicity")); c = _num(gv("toxicity", "max_toxicity_score"))
        st = (PASS if v <= c else FAIL) if (v is not None and c is not None) else NA
        rows.append(_row("Toksisite (ToxinPred2)", v, f"≤ {c}", st,
                         p.methods.get("toxicity", ""), hard=True))

    # ---- IEDB literatür/bilinen-epitop eşleşmesi (pozitif kontrol, bilgilendirici) ----
    ie = m.get("iedb")
    if ie is not None:
        if ie.get("matched") is True:
            org = (ie.get("organisms") or ["?"])[0]
            pmids = ie.get("pmids") or []
            pm = ("PMID " + ", ".join(pmids[:3])) if pmids else "IEDB kaydı"
            val = f"{ie.get('match_type', 'eşleşme')} → {ie.get('epitope_seq', '')} [{org}]"
            rows.append(_row("IEDB bilinen epitop (literatür)", val, pm, PASS,
                             p.methods.get("iedb", ""), kind="doğrulama"))
        elif ie.get("matched") is False:
            rows.append(_row("IEDB bilinen epitop (literatür)", "eşleşme yok",
                             "bilgi (yenilik olabilir)", NA, p.methods.get("iedb", ""),
                             kind="doğrulama"))
        else:
            rows.append(_row("IEDB bilinen epitop (literatür)",
                             ie.get("note", "çalıştırılamadı"), "—", NA,
                             p.methods.get("iedb", ""), kind="doğrulama"))

    return rows
