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


def _pos(location) -> str:
    if not location:
        return "—"
    s = str(location).replace("[", "").replace("]", "")
    return (s.split(":")[0][:9] or "—")


def type_tables(peptides, meta: dict, top_n: int = 15) -> dict:
    """Rapor için TİP BAŞINA (MHC-I / MHC-II / B) sıralı aday tabloları.

    Literatür konvansiyonu: epitoplar tipe göre ayrı tablolarda, her tip kendi
    içinde bağlanma gücüne göre sıralı (T-hücre: %rank artan = güçlü üstte;
    B-hücre: BepiPred azalan). 'star' = tüm zorunlu ölçütleri geçen 'final
    seçilen' epitop: antijenik + güçlü bağlanma (B'de BepiPred geçti) [+ HTL'de
    IFN-γ indükleyici]. Alerjen/toksisite zaten sağ-kalım sert filtresinde
    elendiği için tüm adaylar geçer; yine de sütun gösterilir (literatür gibi).
    top_n + literatürde (IEDB) eşleşen HER aday dahil. Tam döküm Excel'de.
    """
    thr = thr_lookup(meta)
    ag_thr = _num(thr.get(("antigenicity_vaxijen", "threshold"), (0.4,))[0]) or 0.4
    tox_max = _num(thr.get(("toxicity", "max_toxicity_score"), (0.6,))[0]) or 0.6

    def antigenic(p):
        a = _num(p.metrics.get("parent_antigenicity"))
        return a is not None and a >= ag_thr

    def strong(p):
        return any("güçlü" in str(n) for n in getattr(p, "notes", []) or [])

    def matched(p):
        return (p.metrics.get("iedb") or {}).get("matched") is True

    tables: dict[str, list] = {}
    for kind in ("MHC-I", "MHC-II", "B"):
        subs = [p for p in peptides if p.kind == kind]
        if kind == "B":
            subs.sort(key=lambda p: -(_num(p.metrics.get("bcell_score")) or 0))
        else:
            subs.sort(key=lambda p: (_num(p.metrics.get("pseudo_rank")) or 1e9))
        chosen, seen = [], set()
        for p in subs[:top_n]:
            chosen.append(p); seen.add(id(p))
        for p in subs:              # literatür-eşleşenler top_n dışında olsa da
            if id(p) not in seen and matched(p):
                chosen.append(p); seen.add(id(p))
        rows = []
        for p in chosen:
            m = p.metrics
            tox = _num(m.get("toxicity"))
            star = antigenic(p) and (p.passed if kind == "B" else strong(p)) and (
                m.get("ifn_gamma_inducer") is True if kind == "MHC-II" else True)
            rows.append({
                "star": bool(star), "epitope": p.seq, "kind": kind,
                "source": (m.get("locus_tag") or m.get("gene") or p.parent or "?"),
                "pos": _pos(m.get("location")), "length": len(p.seq),
                "rank": _num(m.get("pseudo_rank")), "allele": m.get("best_allele") or "",
                "antigenicity": _num(m.get("parent_antigenicity")),
                "allergen_ok": (m.get("allergen") is not True),
                "toxic_ok": (tox is None or tox < tox_max),
                "immunogenicity": _num(m.get("immunogenicity")),
                "processing": _num(m.get("processing_norm")),
                "bcell": _num(m.get("bcell_score")),
                "ifn_ok": (m.get("ifn_gamma_inducer") is True),
                "iedb": matched(p),
            })
        tables[kind] = rows
    return tables


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
    vir = _num(m.get("virulence"))
    if vir is not None:
        vfid = m.get("vf_identity")
        kw = m.get("vf_keyword")
        ev = (f"VFDB %{vfid}" if m.get("vf_hit")
              else (f"anahtar: {kw}" if kw else "DB kanıtı yok"))
        rows.append(_row("Virülans (VFDB — skorlama, kapı değil)", f"{vir} ({ev})",
                         "bilgi (sert filtre DEĞİL)", NA,
                         m.get("method_discovery", ""), kind="skor"))

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
        if m.get("ifn_gamma_total") is not None:
            v = _num(m.get("ifn_gamma_total")); c = _num(gv("ifn_gamma", "threshold"))
            st = (PASS if v > c else FAIL) if (v is not None and c is not None) else NA
            host = m.get("ifn_gamma_host", "")
            rows.append(_row("IFN-γ indükleme (IFNepitope2)",
                             f"{v}" + (f" [{host}]" if host else ""), f"> {c}", st,
                             p.methods.get("ifn_gamma", "")))
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
