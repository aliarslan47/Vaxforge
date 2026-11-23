"""Adım 8 — Rapor + makine-okur veri paketi + HTML panosu.

Üretilenler (outputs/<run>/):
  candidates.csv     — sıralı aday peptitler + metrikler
  top_peptides.fasta — en iyi adaylar
  run.json           — tüm parametreler/eşikler/özetler (tekrarlanabilirlik)
  report.html        — insan-okur panosu
PDF sonraki sürümde (reportlab/weasyprint).
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from . import correlate
from .i18n import method_label, t
from .models import Peptide


def _candidates_df(peptides: list[Peptide], lang: str = "tr") -> pd.DataFrame:
    from .i18n import method_label
    rows = []
    for i, p in enumerate(peptides, 1):
        rows.append({
            "rank": i, "peptide": p.seq, "kind": p.kind,
            "cds_source": p.parent, "gene": p.metrics.get("gene"),
            "locus_tag": p.metrics.get("locus_tag"), "location": p.metrics.get("location"),
            "candidacy": p.candidacy,
            "mhc_score": p.metrics.get("mhc_score"),
            "best_rank": p.metrics.get("pseudo_rank"),
            "best_allele": p.metrics.get("best_allele"),
            "hosts_presented": ";".join(p.metrics.get("hosts_presented", []) or []),
            "host_coverage": p.metrics.get("host_coverage"),
            "bcell_score": p.metrics.get("bcell_score"),
            "immunogenicity": p.metrics.get("immunogenicity"),
            "processing_norm": p.metrics.get("processing_norm"),
            "conservation_percent": p.metrics.get("conservation_percent"),
            "epitope_conservation": p.metrics.get("epitope_conservation"),
            "self_mimicry": p.metrics.get("self_mimicry"),
            "self_mimicry_host": p.metrics.get("self_mimicry_host"),
            # Dizi karmaşıklığı (bkz. complexity.py) — düşük karmaşıklıklı epitop
            # bağlanma eşiğini geçse de konstrükte konmaz; CSV'de gerekçesi durur.
            "low_complexity": p.metrics.get("low_complexity"),
            "low_complexity_reason": p.metrics.get("low_complexity_reason"),
            "seq_entropy": p.metrics.get("seq_entropy"),
            "seq_max_run": p.metrics.get("seq_max_run"),
            "seq_distinct": p.metrics.get("seq_distinct"),
            "anchor_residues": ";".join(f"{k}={v}" for k, v in
                                        (p.metrics.get("anchor_residues") or {}).items()),
            "allele_anchor_motif": ";".join(f"{k}∈{{{','.join(v)}}}" for k, v in
                                            (p.metrics.get("allele_anchor_motif") or {}).items()),
            "anchor_match": p.metrics.get("anchor_match"),
            # Anchor pozisyonları ALLELE ÖZGÜ ölçülür; hangi tür/allele göre
            # bakıldığı olmadan anchor sütunları yorumlanamaz (2026-11-22).
            "anchor_positions": p.metrics.get("anchor_positions"),
            "anchor_basis": p.metrics.get("anchor_basis"),
            "allergen": p.metrics.get("allergen"),
            "toxicity": p.metrics.get("toxicity"),
            "iedb_match": _iedb_cell(p, "type"),
            "iedb_evidence": _iedb_cell(p, "evidence"),
            "iedb_overlap": _iedb_cell(p, "overlap"),
            "iedb_epitope": _iedb_cell(p, "epitope"),
            "iedb_organism": _iedb_cell(p, "organism"),
            "iedb_assay": _iedb_cell(p, "assay"),
            "iedb_id": _iedb_cell(p, "iedb_id"),
            "iedb_pmids": _iedb_cell(p, "pmids"),
            # YÖNTEM SÜTUNU: eskiden YALNIZ `mhc_score` okunuyordu, bu yüzden
            # B-hücre epitoplarında boş (NaN) kalıyordu — yani run.json/CSV
            # okuyan biri B skorunun GERÇEK BepiPred'den mi yoksa klasik
            # Kolaskar-Tongaonkar+Parker yedeğinden mi geldiğini GÖREMİYORDU.
            # BepiPred sarmalayıcısı yer tutucu yolda kalıp sessizce yedeğe
            # düşerken bu boşluk tam olarak bunu gizledi (DURUM.md bulgu 24).
            # Tipine uygun yöntem anahtarı seçilir, yoksa diğeri denenir.
            "method": method_label(
                p.methods.get("mhc_score") or p.methods.get("bcell_score"), lang),
            "passed": p.passed,
        })
    return pd.DataFrame(rows)


def _iedb_cell(p: Peptide, field: str):
    ie = p.metrics.get("iedb") or {}
    if ie.get("matched") is not True:
        return "" if field != "type" else ("yok" if ie.get("matched") is False else "n/a")
    if field == "type":
        return ie.get("match_type", "eşleşme")
    if field == "evidence":
        return ie.get("evidence", "")
    if field == "overlap":
        return ie.get("overlap", "")
    if field == "epitope":
        return ie.get("epitope_seq", "")
    if field == "organism":
        return (ie.get("organisms") or [""])[0]
    if field == "assay":
        return ";".join(ie.get("assay_evidence") or [])
    if field == "iedb_id":
        return ie.get("iedb_id") or ""
    if field == "pmids":
        return ";".join(ie.get("pmids") or [])
    return ""


def write_package(outdir: Path, peptides: list[Peptide],
                  run_meta: dict) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {}

    # CSV
    lang = run_meta.get("lang", "tr")
    df = _candidates_df(peptides, lang)
    p_csv = outdir / "candidates.csv"
    df.to_csv(p_csv, index=False)
    paths["csv"] = p_csv

    # FASTA (top 20)
    # DÜŞÜK KARMAŞIKLIK bu dosyadan DIŞLANIR: bu dosya sentez/sipariş için
    # kullanılıyor, homopolimerik peptit ise aşı adayı değil (TCR tanıması,
    # konak self-benzerliği, sentezde agregasyon — bkz. complexity.py). Ama
    # SESSİZCE silinmez: dosyanın sonuna EXCLUDED etiketiyle eklenir, böylece
    # "neden yok" sorusu dosyanın içinden yanıtlanır.
    p_fa = outdir / "top_peptides.fasta"
    ok = [p for p in peptides if not p.metrics.get("low_complexity")]
    low = [p for p in peptides if p.metrics.get("low_complexity")][:20]
    with p_fa.open("w") as fh:
        for i, p in enumerate(ok[:20], 1):
            g = p.metrics.get("gene") or ""
            lt = p.metrics.get("locus_tag") or ""
            fh.write(f">cand{i}|{p.kind}|score={p.candidacy}|cds={p.parent}"
                     f"{('|gene='+g) if g else ''}{('|locus='+lt) if lt else ''}\n{p.seq}\n")
        for i, p in enumerate(low, 1):
            fh.write(f">excluded{i}|{p.kind}|score={p.candidacy}|cds={p.parent}"
                     f"|EXCLUDED=low_complexity"
                     f"|reason={p.metrics.get('low_complexity_reason','')}\n{p.seq}\n")
    paths["fasta"] = p_fa

    # Excel (tam liste) — Adaylar + aday-başına TÜM araç sonuçları
    try:
        from . import evaluate
        thr = evaluate.thr_lookup(run_meta)
        tool_rows = []
        for i, p in enumerate(peptides, 1):
            for r in evaluate.candidate_rows(p, thr):
                tool_rows.append({
                    "#": i, "peptit": p.seq, "tip": p.kind, "adaylık": p.candidacy,
                    "araç": r["tool"], "sonuç": r["value"], "eşik": r["cutoff"],
                    "durum": r["status"].replace("✅ ", "").replace("❌ ", ""),
                    "sert_filtre": r["hard"], "yöntem": r["method"],
                })
        p_xlsx = outdir / "candidates_full.xlsx"
        with pd.ExcelWriter(p_xlsx, engine="openpyxl") as xw:
            df.to_excel(xw, sheet_name="Adaylar", index=False)
            pd.DataFrame(tool_rows).to_excel(xw, sheet_name="Arac_sonuclari", index=False)
        paths["xlsx"] = p_xlsx
    except Exception as e:
        (outdir / "xlsx_error.txt").write_text(f"Excel üretilemedi: {e}")

    # Peptit dökümü — raporu KOŞUYU TEKRARLAMADAN yeniden üretebilmek için
    # (bkz. peptide_io; scripts/rebuild_report.py bunu okur).
    try:
        from . import peptide_io
        paths["peptides"] = peptide_io.dump(peptides, outdir / peptide_io.FILENAME)
    except Exception as e:
        (outdir / "peptides_error.txt").write_text(f"peptides.json üretilemedi: {e}")

    # JSON (tam koşum)
    run_meta = dict(run_meta)
    run_meta["candidates"] = df.to_dict(orient="records")
    # Korelat uyarıları JSON'a da makine-okunur yazılır: aşağı akış bir araç
    # bu raporu okuyup "koruyucu aday" iddiası kurmaya kalkarsa, hangi
    # iddianın geçersiz olduğunu buradan görebilir.
    _corr_org = (run_meta.get("organism")
                 or correlate.detect_from_input(run_meta.get("input_path"),
                                                extra=str(run_meta.get("input", ""))))
    run_meta["correlate"] = {
        **correlate.summary(_corr_org),
        "caveats": [
            {"id": c.id, "severity": c.severity, "title": c.title(lang),
             "text": c.text(lang), "blocks_claim": c.blocks_claim,
             "citations": list(c.citations)}
            for c in correlate.caveats(_corr_org)
        ],
    }
    p_json = outdir / "run.json"
    p_json.write_text(json.dumps(run_meta, indent=2, ensure_ascii=False, default=str))
    paths["json"] = p_json

    # HTML panosu
    p_html = outdir / "report.html"
    try:
        p_html.write_text(_html(df, run_meta, peptides), encoding="utf-8")
        paths["html"] = p_html
    except Exception as e:
        (outdir / "html_error.txt").write_text(f"HTML üretilemedi: {e}")

    # PDF (yayın-tarzı) — reportlab varsa
    try:
        from . import report_pdf
        paths["pdf"] = report_pdf.build(outdir, peptides, run_meta)
    except Exception as e:
        (outdir / "pdf_error.txt").write_text(f"PDF üretilemedi: {e}")
    return paths


def _type_tables_html(peptides, meta) -> str:
    """Tip başına (CTL/HTL/B) sıralı, renkli aday tabloları (literatür stili).
    ★ + yeşil satır = tüm zorunlu ölçütleri geçen final seçilen epitop.
    Tam araç dökümü candidates_full.xlsx'te."""
    from . import evaluate
    lang = meta.get("lang", "tr")
    tt = evaluate.type_tables(peptides, meta, top_n=15)
    total = {k: sum(1 for p in peptides if p.kind == k) for k in ("MHC-I", "MHC-II", "B")}
    # Emoji kare yerine CSS rozet (her sistemde düzgün görünür, 'kutu' çıkmaz)
    def _badge(color):
        return (f"<span style='display:inline-block;width:11px;height:11px;border-radius:3px;"
                f"background:{color};vertical-align:middle;margin-right:6px'></span>")
    titles = {
        "MHC-I": f"{_badge('#1565c0')}{t(lang,'tt_ctl')} <span class='mono'>({t(lang,'tt_rank_asc')})</span>",
        "MHC-II": f"{_badge('#b26a00')}{t(lang,'tt_htl')} <span class='mono'>({t(lang,'tt_rank_asc')})</span>",
        "B": f"{_badge('#2e7d32')}{t(lang,'tt_bcell')} <span class='mono'>({t(lang,'tt_bepi_desc')})</span>"}

    def ok(b):
        """ÜÇ DEĞERLİ: None = ölçülmedi. '✗' basmak ölçüm yokluğunu olumsuz
        bulguya çevirirdi (bkz. evaluate.type_tables IFN-γ notu)."""
        if b is None:
            return f'<span class="warn">{t(lang,"ifn_unmeasured")}</span>'
        return '<span class="pass">✔</span>' if b else '<span class="fail">✗</span>'

    def hazard(state, lang_=lang, value=None):
        """Tehlike sütunu: DURUMUN KENDİSİNİ yazar, 'ok' bayrağı basmaz.

        Eskiden burada ✔/✗ vardı ve ✔ "alerjen DEĞİL" demekti; sütun başlığı
        tehlikenin adını taşıdığı için tam tersi okunuyordu. Artık hücre
        'Temiz' / 'Alerjen' / 'Toksik' / '? alan dışı' yazar, toksisitede ham
        skor da parantezde görünür (çıplak işaret bu araç için bilgi taşımıyor).
        """
        txt = t(lang_, f"haz_{state}")
        cls = {"clean": "pass", "allergen": "fail", "toxic": "fail"}.get(state, "warn")
        v = "" if value is None else f' <span class="mono">({value:.2f})</span>'
        return f'<span class="{cls}">{txt}</span>{v}'

    def nz(v, nd=2):
        return "—" if v is None else f"{v:.{nd}f}"

    S = lambda k: t(lang, k)  # noqa: E731
    intro = f'<p style="font-size:.85rem">{t(lang,"tt_intro")}</p>'
    out = [intro]
    for kind in ("MHC-I", "MHC-II", "B"):
        rows = tt.get(kind, [])
        if not rows:
            continue
        out.append(f'<h3>{titles[kind]} <span class="mono">— {len(rows)}/{total[kind]} {S("tt_showing")}</span></h3>')
        if kind == "B":
            head = (f"<tr><th>★</th><th>{S('col_epitope')}</th><th>{S('col_source')}</th><th>{S('col_pos')}</th><th>{S('col_len')}</th>"
                    f"<th>BepiPred</th><th>{S('col_antig')}</th><th>{S('col_allergen')}</th><th>{S('col_toxic')}</th><th>📖</th></tr>")
        elif kind == "MHC-I":
            head = (f"<tr><th>★</th><th>{S('col_epitope')}</th><th>{S('col_source')}</th><th>{S('col_pos')}</th><th>%rank</th>"
                    f"<th>{S('col_allele')}</th><th>{S('col_antig')}</th><th>{S('col_allergen')}</th><th>{S('col_toxic')}</th><th>{S('col_immuno')}</th>"
                    f"<th>{S('col_proc')}</th><th>{S('col_cons')}</th><th>📖</th></tr>")
        else:
            head = (f"<tr><th>★</th><th>{S('col_epitope')}</th><th>{S('col_source')}</th><th>{S('col_pos')}</th><th>%rank</th>"
                    f"<th>{S('col_allele')}</th><th>{S('col_antig')}</th><th>{S('col_allergen')}</th><th>{S('col_toxic')}</th><th>IFN-γ</th>"
                    f"<th>{S('col_cons')}</th><th>📖</th></tr>")
        trs = ""
        for r in rows:
            cls = ' class="sel"' if r["star"] else ""
            star = "★" if r["star"] else ""
            lit = "📖" if r["iedb"] else "–"
            # Düşük karmaşıklık: epitop satırda KALIR, yanında gerekçesi yazar.
            lc = ""
            if r.get("low_complexity"):
                lc = (f'<br><span class="warn" style="font-size:.78rem">⚠ '
                      f'{S("lc_excluded")} — {r.get("low_complexity_reason") or ""}</span>')
            base = (f'<td>{star}</td><td class="mono"><b>{r["epitope"]}</b>{lc}</td>'
                    f'<td class="mono">{r["source"]}</td><td>{r["pos"]}</td>')
            haz = (f'<td>{hazard(r["allergen_state"])}</td>'
                   f'<td>{hazard(r["toxic_state"], value=r.get("toxicity"))}</td>')
            if kind == "B":
                trs += (f'<tr{cls}>{base}<td>{r["length"]}</td><td>{nz(r["bcell"])}</td>'
                        f'<td>{nz(r["antigenicity"])}</td>{haz}<td>{lit}</td></tr>')
            elif kind == "MHC-I":
                trs += (f'<tr{cls}>{base}<td>{nz(r["rank"])}</td><td class="mono">{r["allele"]}</td>'
                        f'<td>{nz(r["antigenicity"])}</td>{haz}<td>{nz(r["immunogenicity"])}</td>'
                        f'<td>{nz(r["processing"])}</td><td>–</td><td>{lit}</td></tr>')
            else:
                trs += (f'<tr{cls}>{base}<td>{nz(r["rank"])}</td><td class="mono">{r["allele"]}</td>'
                        f'<td>{nz(r["antigenicity"])}</td>{haz}<td>{ok(r["ifn_ok"])}</td>'
                        f'<td>–</td><td>{lit}</td></tr>')
        out.append(f"<table>{head}{trs}</table>")
        if kind == "MHC-II" and any(r["ifn_ok"] is None for r in rows):
            out.append(f'<p class="warn" style="font-size:.82rem">⚠ {t(lang,"ifn_note")}</p>')
    return "".join(out)


def _iedb_section(meta: dict, peptides) -> str:
    """IEDB literatür/bilinen-epitop eşleşmesi + (varsa) recall benchmark bölümü."""
    lang = meta.get("lang", "tr")
    im = meta.get("iedb_match")
    if not im:
        return ""
    if not im.get("available"):
        return (f'<h2>{t(lang,"iedb_title")}</h2>'
                f'<div class="warn">{im.get("note", t(lang,"iedb_unavail"))}</div>')
    peptides = peptides or []
    # eşleşen adaylar — doğrulanmış (exact/içerme) önce, sonra zayıf benzerlik
    matched = [(i, p, p.metrics.get("iedb") or {})
               for i, p in enumerate(peptides, 1)
               if (p.metrics.get("iedb") or {}).get("matched") is True]
    matched.sort(key=lambda t3: (0 if t3[2].get("evidence") == "validated" else 1,
                                 -(t3[2].get("overlap") or 0)))
    trs = ""
    for i, p, ie in matched:
        org = "; ".join(o.split(" (")[0] for o in (ie.get("organisms") or [])[:2]) or "—"
        ag = "; ".join((ie.get("antigens") or [])[:2]) or "—"
        assay = "; ".join((ie.get("assay_evidence") or [])[:2]) or "—"
        pmids = ie.get("pmids") or []
        refs = [f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pm}/">{pm}</a>' for pm in pmids[:5]]
        if ie.get("iedb_id"):
            refs.append(f'<a href="https://www.iedb.org/epitope/{ie.get("iedb_id")}">'
                        f'IEDB {ie.get("iedb_id")}</a>')
        plinks = ", ".join(refs) or "—"
        ev = ie.get("evidence")
        strength = t(lang, "strength_validated") if ev == "validated" else t(lang, "strength_weak")
        color = "#1a7f37" if ev == "validated" else "#b26a00"
        trs += (f"<tr><td>#{i}</td><td class='mono'>{p.seq}</td><td>{p.kind}</td>"
                f"<td style='color:{color};font-weight:600'>{strength}</td><td>{ie.get('overlap','')}</td>"
                f"<td>{ie.get('match_type','')}</td><td class='mono'>{ie.get('epitope_seq','')}</td>"
                f"<td>{org}</td><td>{ag}</td><td>{assay}</td><td>{plinks}</td></tr>")
    if not trs:
        trs = f"<tr><td colspan='11'><i>{t(lang,'iedb_none')}</i></td></tr>"
    bm = im.get("benchmark") or {}
    bm_html = ""
    if bm and bm.get("recall") is not None:
        k = bm.get("k", 8)
        dom = ""
        if bm.get("top_ref_share") and bm.get("top_ref_share") >= 0.5:
            dom = (f"<p class='warn' style='font-size:.82rem'>"
                   f"{t(lang,'iedb_ref_dominance').format(share=round(bm['top_ref_share']*100), pmid=bm.get('top_ref') or '—')}</p>")
        rex = bm.get("recall_exact")
        rex_cell = f"<td><b>{round(rex*100,1)}%</b></td>" if rex is not None else "<td>—</td>"
        bm_html = (
            f"<h3>{t(lang,'iedb_val_title')}</h3>"
            f"<p style='font-size:.85rem'>{t(lang,'iedb_val_text').format(k=k)}</p>"
            f"<table><tr><th>{t(lang,'col_known')}</th><th>{t(lang,'col_captured')}</th>"
            f"<th>{t(lang,'col_recall')}</th><th>{t(lang,'col_recall_exact')}</th>"
            f"<th>{t(lang,'col_pred')}</th><th>{t(lang,'col_pred_match')}</th>"
            f"<th>{t(lang,'col_matchrate')}</th></tr>"
            f"<tr><td>{bm['n_known']}</td><td>{bm['n_known_hit']}</td>"
            f"<td><b>{round(bm['recall']*100,1)}%</b></td>{rex_cell}<td>{bm['n_pred']}</td>"
            f"<td>{bm['n_pred_matched']}</td>"
            f"<td><b>{round((bm['precision_like'] or 0)*100,1)}%</b></td></tr></table>{dom}"
        )
    return (
        f"<h2>{t(lang,'iedb_title')}</h2>"
        f"<p style='font-size:.85rem'>{t(lang,'iedb_intro').format(src=im.get('source',''), n_val=im.get('n_validated',0), n_weak=im.get('n_weak',0), tot=im.get('n_candidates','?'))}</p>"
        f"<p class='warn' style='font-size:.82rem'>{t(lang,'iedb_caveat')}</p>"
        f"<table><tr><th>{t(lang,'col_order')}</th><th>{t(lang,'col_candidate')}</th><th>{t(lang,'col_type')}</th>"
        f"<th>{t(lang,'col_strength')}</th><th>{t(lang,'col_overlap')}</th><th>{t(lang,'col_match')}</th>"
        f"<th>{t(lang,'col_iedb_epi')}</th><th>{t(lang,'col_organism')}</th><th>{t(lang,'col_antigen')}</th>"
        f"<th>{t(lang,'col_assay')}</th><th>{t(lang,'col_reference')}</th></tr>"
        f"{trs}</table>{bm_html}"
    )


def _mev_section(meta: dict) -> str:
    """Çok-epitoplu aşı konstruktu (MEV): dizi + özellik tablosu (atıflı)."""
    mv = meta.get("mev")
    if not mv:
        return ""
    lang = meta.get("lang", "tr")
    p = mv["properties"]
    ag, al, tx = p["antigenicity"], p["allergen"], p["toxicity"]
    hs, ss, ds, so = (p["human_similarity"], p.get("secondary_structure", {}),
                      p.get("disorder", {}), p.get("solubility", {}))
    ttl = "Vaccine Construct (MEV)" if lang == "en" else "Aşı Konstruktu (MEV)"
    is_en = lang == "en"

    def yn(v, yes, no):
        return "—" if v is None else (yes if v else no)

    rows = [
        ("Molecular weight" if is_en else "Molekül ağırlığı", f"{p['mw_kda']} kDa"),
        ("Theoretical pI" if is_en else "Teorik pI", p["pI"]),
        ("Negatively charged (Asp+Glu)" if is_en else "Negatif yüklü (Asp+Glu)", p["neg_residues"]),
        ("Positively charged (Arg+Lys)" if is_en else "Pozitif yüklü (Arg+Lys)", p["pos_residues"]),
        ("Instability index", f"{p['instability']} ({'Stable' if p['stable'] else 'Unstable'})"
         if is_en else f"{p['instability']} ({'Kararlı' if p['stable'] else 'Kararsız'})"),
        ("Aliphatic index", p["aliphatic_index"]),
        ("GRAVY", p["gravy"]),
        ("Antigenicity", f"{ag['score']} ({yn(ag['antigenic'],'Antigenic' if is_en else 'Antijenik','Non' )}) · {ag['method']}"),
        ("Allergenicity", yn(al["allergenic"], "Allergen", "Non-allergen") + f" · {al['method']}"
         + (f" · {al['n_matches']} 6-mer" if al.get("n_matches") else "")),
        ("Toxicity", yn(tx["toxic"], "Toxic", "Non-toxic") + f" · {tx['method']}"),
        ("Human proteome similarity" if is_en else "İnsan proteomu benzerliği",
         ("No similarity" if is_en else "Benzerlik yok") if hs.get("available") and not hs.get("similar")
         else (f"{hs.get('best_pident')}%" if hs.get("available") else "—")),
        ("Solubility" if is_en else "Çözünürlük",
         f"{so.get('scaled_solubility')} ({yn(so.get('soluble'),'Soluble' if is_en else 'Çözünür','Low')})"
         if so.get("scaled_solubility") is not None else "—"),
        ("Disordered regions" if is_en else "Düzensiz bölgeler",
         f"{ds.get('percent_disordered')}%" if ds.get("percent_disordered") is not None else "—"),
        ("α-Helix", f"{ss.get('helix_pct')}%" if ss.get("helix_pct") is not None else "—"),
        ("β-Strand", f"{ss.get('strand_pct')}%" if ss.get("strand_pct") is not None else "—"),
        ("Coil", f"{ss.get('coil_pct')}%" if ss.get("coil_pct") is not None else "—"),
    ]
    trows = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    nb = mv["n_by_kind"]
    comp = (f"β-defensin adjuvant + EAAAK/AAY/GPGPG/KK linkers · "
            f"CTL={nb['MHC-I']}, HTL={nb['MHC-II']}, LBL={nb['B']}" if is_en
            else f"β-defensin adjuvan + EAAAK/AAY/GPGPG/KK linker · "
                 f"CTL={nb['MHC-I']}, HTL={nb['MHC-II']}, LBL={nb['B']}")
    seq_lbl = "Construct sequence" if is_en else "Konstrukt dizisi"
    len_lbl = "length" if is_en else "uzunluk"
    # Düşük karmaşıklık nedeniyle kasete KONMAYAN epitoplar: bağlanma eşiğini
    # geçtikleri için "neden yok" sorusu doğar, o yüzden gerekçesiyle listelenir.
    exc = mv.get("excluded_low_complexity") or []
    exc_html = ""
    if exc:
        erows = "".join(
            f'<tr><td class="mono">{e["seq"]}</td><td>{e["kind"]}</td>'
            f'<td class="mono">{e.get("source","")}</td>'
            f'<td>{"—" if e.get("rank") is None else e["rank"]}</td>'
            f'<td>{e.get("reason","")}</td></tr>' for e in exc)
        exc_html = (
            f"<h3>⚠ {t(lang,'lc_title')} <span class='mono'>— {len(exc)}</span></h3>"
            f"<p style='font-size:.85rem'>{t(lang,'lc_intro')}</p>"
            f"<table><tr><th>{t(lang,'col_epitope')}</th><th>{'Type' if is_en else 'Tip'}</th>"
            f"<th>{t(lang,'col_source')}</th><th>%rank</th>"
            f"<th>{'Reason' if is_en else 'Gerekçe'}</th></tr>{erows}</table>")
    # KONSTRÜKT ALERJENİTESİ — uyarı kutusu. Adaylarda sert filtre olduğu için
    # her epitop temizdir; konstrükt yine alerjen çıkıyorsa eşleşme birleşme
    # yerinde doğmuş olabilir ve LİNKER DEĞİŞTİRİLEREK giderilebilir. Eskiden
    # bu satır tabloda tek kelimeydi ("Allergen") ve eyleme dönüşmüyordu.
    al_html = ""
    if al.get("allergenic") and al.get("matches"):
        reg_lbl = {"junction": ("birleşme yeri (epitop/linker)", "junction (epitope/linker)"),
                   "linker": ("linker dizisi", "linker sequence"),
                   "adjuvant": ("adjuvan (β-defensin, konak kökenli — beklenir)",
                                "adjuvant (β-defensin, host-derived — expected)"),
                   "epitope": ("EPİTOBUN KENDİSİ — aday filtresi kaçırmış",
                               "THE EPITOPE ITSELF — candidate filter missed it"),
                   "bilinmiyor": ("bilinmiyor", "unknown")}
        arows = "".join(
            f'<tr><td class="mono">{h["kmer"]}</td><td>{h["start"] + 1}</td>'
            f'<td>{reg_lbl.get(h["region"], (h["region"], h["region"]))[1 if is_en else 0]}</td>'
            f'<td class="mono">{h["detail"]}</td>'
            f'<td class="mono">{h.get("source", "")}</td></tr>' for h in al["matches"])
        if al.get("from_epitope"):
            verdict = ("At least one match lies INSIDE an epitope: that epitope should be "
                       "dropped from the construct (the candidate-level filter missed it)."
                       if is_en else
                       "En az bir eşleşme bir EPİTOBUN İÇİNDE: o epitop konstrüktten "
                       "çıkarılmalı (aday düzeyi filtre kaçırmış).")
            vcls = "fail"
        elif al.get("fixable_by_linker"):
            verdict = ("All matches arise at linkers or at epitope/linker junctions — the "
                       "epitopes themselves are clean. This is fixable by changing the "
                       "linker or the epitope order; it is not a property of the antigens."
                       if is_en else
                       "Eşleşmelerin tamamı linkerlarda ya da epitop/linker birleşme "
                       "yerlerinde doğuyor — epitopların kendisi temiz. Linker ya da epitop "
                       "sırası değiştirilerek giderilebilir; antijenlerin özelliği değil.")
            vcls = "warn"
        else:
            verdict = ("Matches come from the adjuvant, which is host-derived; expected, "
                       "not a filtering criterion." if is_en else
                       "Eşleşmeler konak kökenli adjuvandan geliyor; beklenen durum, "
                       "eleme ölçütü değil.")
            vcls = "warn"
        al_html = (
            f"<h3>⚠ {'Construct allergenicity (FAO/WHO 6-mer)' if is_en else 'Konstrükt alerjenitesi (FAO/WHO 6-mer)'}"
            f" <span class='mono'>— {al['n_matches']}</span></h3>"
            f"<p style='font-size:.85rem'><span class='{vcls}'>{verdict}</span></p>"
            f"<table><tr><th>6-mer</th><th>{'Position' if is_en else 'Konum'}</th>"
            f"<th>{'Origin' if is_en else 'Köken'}</th><th>{'Parts' if is_en else 'Parçalar'}</th>"
            f"<th>{'Matched allergen' if is_en else 'Eşleşen alerjen'}</th></tr>{arows}</table>")
    return (
        f"<h2>🧩 {ttl}</h2>"
        f"<p style='font-size:.85rem'>{comp} · {len_lbl}: <b>{len(mv['sequence'])} aa</b> "
        f"(adjuvan: {mv['adjuvant']})</p>"
        f"<details><summary>{seq_lbl}</summary>"
        f"<p class='mono'>{mv['sequence']}</p></details>"
        f"<table><tr><th>{t(lang,'col_param') if False else ('Property' if is_en else 'Özellik')}</th>"
        f"<th>{'Value' if is_en else 'Değer'}</th></tr>{trows}</table>"
        f"{al_html}{exc_html}"
    )


def _immuno_section(meta: dict) -> str:
    """İmmünojenisite — in silico eleme (ImmForge). DÜRÜST etiketli (aşırı-iddia yok)."""
    imm = meta.get("immunogenicity")
    if not imm or not imm.get("per_host"):
        return ""
    is_en = meta.get("lang", "tr") == "en"
    ttl = "Immunogenicity — in silico screen (ImmForge)" if is_en else "İmmünojenisite — in silico eleme (ImmForge)"
    disc = ("Mechanistic engine PREDICTION — not in vivo efficacy, NOT validated. Only the kinetic shape "
            "was validated on peptide-vaccine data; verdict/isotype/profile are operational (threshold "
            "WORKS≥60). IgM suppressed by single-step class switching." if is_en
            else "Mekanistik motor ÖNGÖRÜSÜ — in vivo etkinlik DEĞİL, valide EDİLMEMİŞ. Yalnız kinetik şekil "
            "peptid-aşısı verisinde doğrulandı; verdict/profil operasyonel (eşik OLUR≥60). İzotip KİMLİKLERİ "
            "literatür-steering'den (IgM primer→IgG alt-sınıf karışımı); RELATİF büyüklükler operasyonel "
            "(Th1/Th2 dengesine bağlı, per-konstrukt valide değil).")
    vmap = {"OLUR": ("WORKS" if is_en else "OLUR"), "ZAYIF": ("WEAK" if is_en else "ZAYIF"),
            "OLMAZ": ("FAILS" if is_en else "OLMAZ")}
    hdr = (["Host", "Verdict", "Score", "MHC-I", "MHC-II", "Clearance", "Profile", "Ig", "Calibration"] if is_en
           else ["Konak", "Karar", "Skor", "MHC-I", "MHC-II", "Temizlik", "Profil", "Ig", "Kalibrasyon"])
    trs = ""
    for ph in imm["per_host"]:
        pr = ph.get("profile") or {}
        prof = f"{pr.get('polarization', '—')} ({round(pr.get('th1', 0) * 100)}/{round(pr.get('th2', 0) * 100)})"
        cal_sp = ph.get("calibration_species", "?")
        cal = (cal_sp if ph.get("species_calibrated")
               else (f"{cal_sp}-proxy" if is_en else f"{cal_sp}-proxy"))
        trs += ("<tr>"
                f"<td>{ph.get('host_label', ph.get('host', ''))}</td>"
                f"<td>{vmap.get(ph.get('verdict', ''), ph.get('verdict', ''))}</td>"
                f"<td>{round(ph.get('score', 0))}/100</td>"
                f"<td>{ph.get('n_mhci', '')}</td><td>{ph.get('n_mhcii', '')}</td>"
                f"<td>{round(ph.get('cleared_fraction', 0) * 100)}%</td>"
                f"<td>{prof}</td><td>{pr.get('dominant_isotype') or '—'}</td><td>{cal}</td></tr>")
    r0 = imm["per_host"][0].get("reasons") or []
    reasons_html = ("<p style='font-size:.8rem;color:#555'>· " + " · ".join(r0) + "</p>") if r0 else ""
    return (
        f"<h2>🛡️ {ttl}</h2>"
        f"<p style='font-size:.8rem;background:#fff7e6;border:1px solid #f0c36d;padding:.5rem;border-radius:6px'>⚠ {disc}</p>"
        "<table><tr>" + "".join(f"<th>{h}</th>" for h in hdr) + f"</tr>{trs}</table>{reasons_html}"
    )


def _structure_section(meta: dict) -> str:
    """Yapısal doğrulama — ESMFold (pLDDT) + RSA + konformasyonel yüzey B (Faz 3)."""
    st = meta.get("structure")
    if not st:
        return ""
    is_en = meta.get("lang", "tr") == "en"
    ttl = "Structural validation — ESMFold" if is_en else "Yapısal doğrulama — ESMFold"
    if not st.get("available"):
        note = st.get("note", "")
        body = (f"<p><i>{note}</i></p>" if note else
                ("<p><i>Structure not predicted.</i></p>" if is_en else "<p><i>Yapı öngörülmedi.</i></p>"))
        return f"<h2>🧬 {ttl}</h2>{body}"
    surf = st.get("surface") or {}
    cb = st.get("conformational_b") or {}
    disc = ("ESMFold PREDICTION via ESMAtlas API — not an experimental structure. "
            "pLDDT>70 confident, <50 disordered. Conformational-B is DiscoTope-LIKE "
            "(RSA+Parker), not the official DiscoTope." if is_en else
            "ESMAtlas API üzerinden ESMFold ÖNGÖRÜSÜ — deneysel yapı değil. pLDDT>70 "
            "güvenilir, <50 düzensiz. Konformasyonel-B, DiscoTope-BENZERİ (RSA+Parker); "
            "resmi DiscoTope değil.")
    # pLDDT GÜVEN BANNER'I: eşik altındaysa yapıdan TÜREYEN sonuçlar (RSA/yüzey,
    # konformasyonel B) "düşük güven" etiketi taşır. MEV linkerlarla birleşmiş
    # yapay bir kimera olduğundan düşük pLDDT BEKLENEN durumdur — kusur değil;
    # banner bunu da yazar, yoksa sayı "kötü sonuç" gibi okunuyor.
    cl = st.get("plddt_confidence") or {}
    thr = st.get("plddt_threshold", 70.0)
    banner = ""
    if cl and not cl.get("reliable", True):
        lbl = cl.get("label_en" if is_en else "label_tr", "")
        banner = (
            f"<p style='font-size:.85rem'><span class='warn'><b>⚠ pLDDT "
            f"{st.get('plddt_mean','—')} &lt; {thr} — {lbl}.</b></span> "
            + ("Everything read OFF this model is therefore low-confidence too: "
               "surface % (RSA) and the conformational B-cell patches below. "
               "Note this is EXPECTED for a MEV: the construct is a linker-joined "
               "chimera, not a natural folding unit, so a low pLDDT is not by itself "
               "a defect — but it means the surface/conformational results cannot "
               "carry the weight of a measurement." if is_en else
               "Dolayısıyla bu modelden OKUNAN her şey de düşük güvenlidir: "
               "aşağıdaki yüzey % (RSA) ve konformasyonel B-hücre yamaları. "
               "Bunun MEV için BEKLENEN olduğunu da not edelim: konstrükt "
               "linkerlarla birleştirilmiş yapay bir kimera, doğal bir katlanma "
               "birimi değil; yani düşük pLDDT tek başına kusur değil — ama "
               "yüzey/konformasyonel sonuçların ölçüm ağırlığı taşımadığı anlamına "
               "geliyor.") + "</p>")
    lc_tag = ("" if cl.get("reliable", True) else
              (" <span class='warn mono'>(low confidence)</span>" if is_en
               else " <span class='warn mono'>(düşük güven)</span>"))
    rows = [
        (("Source" if is_en else "Kaynak"), st.get("source", "—")),
        (("Residues" if is_en else "Kalıntı"), st.get("n_res", "—")),
        ("pLDDT (mean/min)",
         f"{st.get('plddt_mean','—')} / {st.get('plddt_min','—')}"
         + (f" — {cl.get('label_en' if is_en else 'label_tr','')}" if cl else "")),
        (("Surface %" if is_en else "Yüzey %"),
         f"{surf.get('percent_surface','—')}%{lc_tag}"),
        (("Conformational B patches" if is_en else "Konformasyonel B yaması"),
         f"{cb.get('n_patches', '—')}{lc_tag}"),
    ]
    trs = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    patch_html = ""
    for p in (cb.get("patches") or [])[:6]:
        # Her yamanın KENDİ pLDDT'si: güven dizi boyunca eşit dağılmaz, bir yama
        # güvenilir bölgede, başkası düzensiz bölgede olabilir.
        _pl = p.get("plddt_mean")
        _flag = (" <span class='warn'>⚠</span>" if p.get("low_confidence") else "")
        patch_html += (f"<tr><td>{p['start']}–{p['end']}</td><td>{p['length']}</td>"
                       f"<td class='mono'>{p['peptide']}</td><td>{p['mean_score']}</td>"
                       f"<td>{'—' if _pl is None else _pl}{_flag}</td></tr>")
    if patch_html:
        ph = ("<p style='margin-top:.6rem'><b>Surface B-cell patches</b></p>" if is_en
              else "<p style='margin-top:.6rem'><b>Yüzey B-hücre yamaları</b></p>")
        hh = (["Range", "Len", "Peptide", "Score", "pLDDT"] if is_en
              else ["Aralık", "Uzunluk", "Peptit", "Skor", "pLDDT"])
        patch_html = ph + "<table><tr>" + "".join(f"<th>{h}</th>" for h in hh) + f"</tr>{patch_html}</table>"
    return (
        f"<h2>🧬 {ttl}</h2>"
        f"<p style='font-size:.8rem;background:#eef6ff;border:1px solid #9cc3e6;padding:.5rem;border-radius:6px'>ℹ {disc}</p>"
        f"{banner}<table>{trs}</table>{patch_html}"
    )


def _glyco_section(meta: dict) -> str:
    """Glikozilasyon — glikan-maskeli B-epitop filtresi (Faz 4)."""
    g = meta.get("glycosylation")
    if not g or not g.get("available"):
        return ""
    is_en = meta.get("lang", "tr") == "en"
    ttl = "Glycosylation — glycan-masked epitopes" if is_en else "Glikozilasyon — glikan-maskeli epitoplar"
    disc = ("N-X-S/T is an exact sequon motif, but actual glycan occupancy needs "
            "experimental proof; O-glyco here is a Ser/Thr-rich heuristic. PREDICTION — "
            "masked B-cell epitopes may be antibody-inaccessible in vivo (soft flag, not filtered)."
            if is_en else
            "N-X-S/T kesin bir sequon motifidir; gerçek glikan işgali deneysel doğrulanmalı. "
            "O-gliko burada Ser/Thr-yoğun heuristiktir. ÖNGÖRÜ — maskeli B-epitopları in vivo "
            "antikora erişilemeyebilir (yumuşak uyarı; aday elenmedi).")
    rows = [
        (("Candidates scanned" if is_en else "Taranan aday"), g.get("n_candidates", "—")),
        (("B-cell epitopes" if is_en else "B-hücre epitobu"), g.get("n_bcell", "—")),
        (("B-cell glycan-masked" if is_en else "B-hücre glikan-maskeli"), g.get("n_bcell_masked", "—")),
        (("All masked (N-sequon)" if is_en else "Tüm maskeli (N-sequon)"), g.get("n_masked", "—")),
    ]
    trs = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    return (
        f"<h2>🍬 {ttl}</h2>"
        f"<p style='font-size:.8rem;background:#fff7e6;border:1px solid #f0c36d;padding:.5rem;border-radius:6px'>⚠ {disc}</p>"
        f"<table>{trs}</table>"
    )


def _conservation_section(meta: dict) -> str:
    """Suşlar-arası konservasyon — otomatik çok-suş dahil (Faz 8)."""
    c = meta.get("conservation")
    if not c:
        return ""
    is_en = meta.get("lang", "tr") == "en"
    ttl = "Cross-strain conservation" if is_en else "Suşlar-arası konservasyon"
    if not c.get("computed"):
        auto = meta.get("conservation_auto_strains") or {}
        note = auto.get("note") or ("no strain data" if is_en else "suş verisi yok")
        body = (f"<p><i>{'Not computed' if is_en else 'Hesaplanmadı'} — {note}.</i></p>")
        return f"<h2>🧬 {ttl}</h2>{body}"
    rows = [
        (("Strains" if is_en else "Suş sayısı"), c.get("n_strains", "—")),
        (("Conserved candidates" if is_en else "Korunmuş aday"), c.get("n_conserved", "—")),
        (("Threshold" if is_en else "Eşik"), f"≥%{c.get('min_percent', '—')}"),
        # 2026-11-23: bu üç alan csum içinde VARDI ama hiçbir raporda
        # görünmüyordu — yani "kaç adayda ölçüm YAPILAMADI" gizli kalıyordu.
        (("Not measurable (fewer than min. orthologs)" if is_en
          else "Ölçülemedi (en az ortolog sayısının altında)"),
         c.get("n_unmeasured", "—")),
        (("Min. orthologs required" if is_en else "Gereken en az ortolog"),
         c.get("min_orthologs", "—")),
        (("Strains identical to input (excluded)" if is_en
          else "Girdiyle aynı suş (sayılmadı)"), c.get("n_identical_strains", "—")),
        (("Method" if is_en else "Yöntem"), c.get("method", "—")),
    ]
    if c.get("auto_fetched"):
        rows.append((("Auto-fetched" if is_en else "Otomatik çekilen"), c.get("auto_note", "—")))
    trs = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    return f"<h2>🧬 {ttl}</h2><table>{trs}</table>"


def _mimicry_section(meta: dict) -> str:
    """Moleküler mimikri — tam + BLOSUM62 kısmi konak self-eşleşme (Faz 7)."""
    m = meta.get("mimicry")
    if not m or not m.get("available"):
        return ""
    is_en = meta.get("lang", "tr") == "en"
    ttl = "Molecular mimicry — host self-match" if is_en else "Moleküler mimikri — konak self-eşleşme"
    disc = ("Exact k-mer match = strong autoimmunity risk; BLOSUM62 Hamming-1 (single "
            "conservative substitution) = partial cross-reactivity risk. SOFT — flags, "
            "does not filter." if is_en else
            "Tam k-mer eşleşme = güçlü otoimmünite riski; BLOSUM62 Hamming-1 (tek konservatif "
            "ikame) = kısmi çapraz-reaktivite riski. YUMUŞAK — işaretler, elemez.")
    rows = [
        (("k-mer length" if is_en else "k-mer uzunluğu"), m.get("k", "—")),
        (("Candidates" if is_en else "Aday"), m.get("n_candidates", "—")),
        (("Exact mimics (hard)" if is_en else "Tam mimik (sert)"), m.get("n_exact", "—")),
        (("Partial mimics (BLOSUM Hamming-1)" if is_en else "Kısmi mimik (BLOSUM Hamming-1)"),
         m.get("n_near", "—")),
        (("Hosts" if is_en else "Konaklar"), ", ".join(m.get("hosts", [])) or "—"),
    ]
    trs = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    return (
        f"<h2>🧩 {ttl}</h2>"
        f"<p style='font-size:.8rem;background:#fff7e6;border:1px solid #f0c36d;padding:.5rem;border-radius:6px'>⚠ {disc}</p>"
        f"<table>{trs}</table>"
    )


def _diagnostics_section(meta: dict) -> str:
    """Tanılama — yakalanan istisnaların görünür kaydı (iki dilli).

    Pipeline aşamaları try/except ile sarılı olduğu için bir araç çöktüğünde koşu
    devam eder. Bu bölüm olmadan "araç kurulu değil" ile "araç çöktü" raporda
    AYNI görünür. Beklenmeyen hata varsa uyarı, yoksa dürüst temiz notu basar.
    """
    items = meta.get("diagnostics") or []
    summ = meta.get("diagnostics_summary") or {}
    if not items and not summ:
        return ""
    lang = "en" if meta.get("lang", "tr") == "en" else "tr"
    n_unexpected = int(summ.get("unexpected") or 0)
    ttl = t(lang, "diagnostics")

    if not n_unexpected:
        body = ("<p style='font-size:.8rem;background:#eef9f0;border:1px solid #8fc9a0;"
                f"padding:.45rem;border-radius:6px'>✓ {t(lang, 'diag_clean')}</p>")
    else:
        stages = ", ".join(summ.get("stages") or [])
        rows = "".join(
            f"<tr><td>{d.get('stage','—')}</td><td>{d.get('error','—')}</td>"
            f"<td>{(d.get('message') or '—')[:160]}</td><td>{d.get('where') or '—'}</td></tr>"
            for d in items if not d.get("expected"))
        hdr = (("Stage", "Error", "Message", "Where") if lang == "en"
               else ("Aşama", "Hata", "Mesaj", "Yer"))
        body = (
            "<p style='font-size:.8rem;background:#fff6e5;border:1px solid #e0b25c;"
            f"padding:.45rem;border-radius:6px'>⚠ {t(lang, 'diag_unexpected')}: "
            f"{n_unexpected} ({stages})</p>"
            f"<table><tr><th>{hdr[0]}</th><th>{hdr[1]}</th><th>{hdr[2]}</th>"
            f"<th>{hdr[3]}</th></tr>{rows}</table>")
    return f"<h2>🩺 {ttl}</h2>{body}"


def _tool_status_section(meta: dict) -> str:
    """Araç durumu — hangi dış araç bu koşuda GERÇEKTEN hazırdı (iki dilli).

    NEDEN: Referanslar/Yöntemler tablosu `citations.py`'den koşulsuz besleniyordu.
    Kuduz koşusunda işleme omurgasının binary'si segfault verdi; 85 adayın
    tamamında `processing_norm` None kaldı, ama rapor o aracı atfıyla yöntem
    olarak gösterdi. Bu bölüm eksik aracı açıkça adlandırır ve eksik metriğin
    olumsuz bir bulgu OLMADIĞINI söyler.
    """
    status = meta.get("tool_status") or {}
    summ = meta.get("tool_status_summary") or {}
    if not status:
        return ""
    lang = "en" if meta.get("lang", "tr") == "en" else "tr"
    ttl = t(lang, "tool_status")
    miss = ((summ.get("missing_labels_en") or summ.get("missing_labels") or [])
            if lang == "en" else (summ.get("missing_labels") or []))
    count = (f"{summ.get('available', 0)} / {summ.get('total', 0)} "
             f"({t(lang, 'tool_status_count')})")
    if not miss:
        body = ("<p style='font-size:.8rem;background:#eef9f0;border:1px solid #8fc9a0;"
                f"padding:.45rem;border-radius:6px'>✓ {t(lang, 'tool_status_ok')} "
                f"[{count}]</p>")
    else:
        # EKSİKLİĞİN TÜRÜ: "kurulu değil" kalıcı; "doğrulama koşusu düştü"
        # geçici yük hatası olabilir (bir kez retry edilip yine düştü). İkisini
        # aynı etiketle göstermek, çalışan bir aracı yokmuş gibi raporluyordu.
        sfx = "_en" if lang == "en" else ""
        vf = set(summ.get("verify_failed" + sfx) or [])
        pe = set(summ.get("probe_error" + sfx) or [])
        st_vf = ("file present, verification run FAILED (may be transient load)"
                 if lang == "en" else
                 "dosya var, DOĞRULAMA KOŞUSU DÜŞTÜ (geçici yük olabilir)")
        st_pe = ("probe error" if lang == "en" else "kontrol hatası")
        rows = "".join(
            f"<tr><td>{lbl}</td><td>"
            + (f"<span class='warn'><b>{st_vf}</b></span>" if lbl in vf
               else f"<span class='warn'>{st_pe}</span>" if lbl in pe
               else f"<b>{t(lang, 'tool_status_absent')}</b>")
            + "</td></tr>"
            for lbl in miss)
        body = (
            "<p style='font-size:.8rem;background:#fff6e5;border:1px solid #e0b25c;"
            f"padding:.45rem;border-radius:6px'>⚠ {t(lang, 'tool_status_warn')} "
            f"[{count}]</p>"
            f"<table><tr><th>{t(lang, 'tool_status_tool')}</th>"
            f"<th>{t(lang, 'tool_status_state')}</th></tr>{rows}</table>")
    # Kontrol DERİNLİĞİ de raporlanır: "dosya var" bir garanti değil.
    nv, ns = int(summ.get("verified") or 0), int(summ.get("shallow") or 0)
    if nv or ns:
        body += ("<p style='font-size:.75rem;color:#555;margin-top:.3rem'>"
                 + t(lang, "tool_check_note").format(v=nv, s=ns) + "</p>")
    return f"<h2>🧰 {ttl}</h2>{body}"


def _cloning_section(meta: dict) -> str:
    """In silico klonlama — kodon optimizasyonu + CAI + restriksiyon + vektör (Faz 5)."""
    cl = meta.get("cloning")
    if not cl or not cl.get("available"):
        return ""
    is_en = meta.get("lang", "tr") == "en"
    ttl = "In silico cloning — codon optimization" if is_en else "In silico klonlama — kodon optimizasyonu"
    disc = ("Codon-optimized for E. coli (Kazusa); CAI/GC computed. PREDICTION — not a "
            "guarantee of expression. Internal restriction sites break classic cloning." if is_en else
            "E. coli için kodon-optimize (Kazusa); CAI/GC hesaplandı. ÖNGÖRÜ — ekspresyon "
            "garantisi değil. İç restriksiyon bölgeleri klasik klonlamayı bozar.")
    sites = cl.get("restriction_sites") or {}
    sites_str = ", ".join(f"{e}×{len(p)}" for e, p in sites.items()) or ("none" if is_en else "yok")
    lang = "en" if is_en else "tr"
    rows = [
        (("Host" if is_en else "Konak"), cl.get("host", "—")),
        (("Protein length (aa)" if is_en else "Protein (aa)"), cl.get("protein_len", "—")),
        (("DNA length (bp)" if is_en else "DNA (bp)"), cl.get("dna_len_bp", "—")),
        (t(lang, "cai_optimized"), cl.get("cai", "—")),
        (t(lang, "cai_baseline"), cl.get("cai_baseline", "—")),
        ("GC %", f"{cl.get('gc_percent', '—')}%"),
        (("Internal restriction sites" if is_en else "İç restriksiyon bölgeleri"), sites_str),
        (("Suggested vector" if is_en else "Önerilen vektör"), cl.get("suggested_vector", "—")),
        (("Tag / strategy" if is_en else "Etiket / strateji"), cl.get("tag", "—")),
    ]
    n_skip = cl.get("n_skipped_residues") or 0
    if n_skip:
        rows.append((t(lang, "cai_skipped"), n_skip))
    trs = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    # CAI totoloji uyarısı HER ZAMAN gösterilir (iki dilli); atlanan kalıntı varsa ek uyarı.
    warns = [t(lang, "cai_tautology_warn")]
    if n_skip:
        warns.append(t(lang, "cai_skipped_warn"))
    warn_html = "".join(
        "<p style='font-size:.78rem;background:#fff6e5;border:1px solid #e0b25c;"
        f"padding:.45rem;border-radius:6px;margin:.3rem 0'>⚠ {w}</p>" for w in warns)
    dna = cl.get("dna", "")
    dna_prev = (dna[:120] + "…") if len(dna) > 120 else dna
    return (
        f"<h2>🧫 {ttl}</h2>"
        f"<p style='font-size:.8rem;background:#eef6ff;border:1px solid #9cc3e6;padding:.5rem;border-radius:6px'>ℹ {disc}</p>"
        f"<table>{trs}</table>"
        f"{warn_html}"
        f"<p class='mono' style='margin-top:.4rem'>5'-{dna_prev}-3'</p>"
    )


def _methods_section(meta: dict) -> str:
    """Kullanılan yöntemler — funnel özetinden GERÇEK/vekil araç etiketleri (Faz 1).

    Lokalizasyonun PSORTb (GERÇEK) mı yoksa heuristik mi olduğunu; TM/sinyal/
    antijenite/öz-homoloji araçlarını açıkça gösterir. Bilim değişmez — yalnız
    'hangi araç koştu' şeffaf raporlanır (uydurma yok — çekirdek ilke).
    """
    fm = meta.get("funnel_methods") or {}
    if not fm:
        return ""
    lang = meta.get("lang", "tr")
    rows = [
        ("meth_localization", fm.get("lokalizasyon")),
        ("meth_tm", fm.get("transmembran")),
        ("meth_signal", fm.get("sinyal_peptidi")),
        ("meth_antigenicity", fm.get("antijenite")),
        ("meth_antigenicity_vaxijen", fm.get("antijenite_caprazkontrol")),
        ("meth_homology", fm.get("self_homoloji")),
    ]
    trs = "".join(
        f"<tr><td>{t(lang, key)}</td><td>{method_label(str(val), lang)}</td></tr>"
        for key, val in rows if val
    )
    if not trs:
        return ""
    # Bakteri koşusunda izinli lokalizasyon listesi Gram'a göre seçilir; hangi
    # mantığın uygulandığı rapora yazılır (tekrarlanabilirlik + iki dilli).
    gram_note = ""
    if meta.get("profile") == "bacteria" and meta.get("gram") in ("positive", "negative"):
        gram_note = ("<p style='font-size:.78rem;background:#eef6ff;border:1px solid #9cc3e6;"
                     f"padding:.45rem;border-radius:6px;margin-top:.3rem'>ℹ {t(lang, 'gram_allowed_note')}</p>")
    return (
        f"<h2>🧪 {t(lang, 'rep_methods')}</h2>"
        f"<p style='font-size:.85rem'>{t(lang, 'rep_methods_sub')}</p>"
        f"<table><tr><th>{t(lang, 'meth_step')}</th><th>{t(lang, 'meth_method')}</th></tr>{trs}</table>"
        f"{gram_note}"
    )


def _correlate_section(meta: dict) -> str:
    """İMMÜN KORELAT UYARILARI — organizmaya bağlı, ZORUNLU bölüm.

    Bu pipeline epitop ÖNGÖRÜSÜ üretir; korumayı belirleyen ise patojenin
    KORELATIDIR. Bazı organizmalarda korelat protein tabanlı bir pipeline'ın
    göremeyeceği bir moleküle bağlıdır (Leptospira'da LPS). Bu bölüm her
    organizmada yazılır; hiçbir koşulda sessiz kalınmaz.
    """
    lang = meta.get("lang", "tr")
    org = (meta.get("organism")
           or correlate.detect_from_input(meta.get("input_path"),
                                          extra=str(meta.get("input", ""))))
    cavs = correlate.caveats(org)
    if not cavs:
        return ""
    colors = {"critical": "#b00020", "major": "#9a6700", "info": "#555"}
    bg = {"critical": "#fdecef", "major": "#fff8e6", "info": "#f5f5f5"}
    blocks = []
    if any(c.severity == "critical" for c in cavs):
        blocks.append(f'<div class="warn"><b>{t(lang, "corr_critical")}</b></div>')
    for c in cavs:
        cits = "".join(f"<li>{x}</li>" for x in c.citations)
        claim = (f'<p><i>{t(lang, "corr_blocks")} &ldquo;{c.blocks_claim}&rdquo;</i></p>'
                 if c.blocks_claim else "")
        blocks.append(
            f'<div style="border-left:4px solid {colors.get(c.severity, "#555")};'
            f'background:{bg.get(c.severity, "#f5f5f5")};padding:.5rem .8rem;'
            f'margin:.5rem 0;border-radius:4px">'
            f'<p style="margin:.1rem 0"><b style="color:{colors.get(c.severity, "#555")}">'
            f'[{c.severity.upper()}]</b> <b>{c.title(lang)}</b></p>'
            f'<p>{c.text(lang)}</p>{claim}'
            + (f'<ul class="mono">{cits}</ul>' if cits else "")
            + '</div>')
    return (f'<h2>{t(lang, "corr_title")}</h2>'
            f'<p>{t(lang, "corr_intro")}</p>' + "".join(blocks))


def _curated_validation_section(meta: dict) -> str:
    """Organizmaya özgü küratörlü doğrulama seti (IEDB boş olduğunda)."""
    cv = meta.get("curated_validation") or {}
    if not cv:
        return ""
    lang = meta.get("lang", "tr")
    tr = lang != "en"
    hc = cv.get("host_coverage") or {}
    host_rows = "".join(
        f"<tr><td>{h}</td><td>{d.get('n_total', 0)}</td>"
        f"<td>{d.get('n_positive', 0)}</td><td>{d.get('n_negative', 0)}</td></tr>"
        for h, d in sorted(hc.items(), key=lambda kv: -kv[1].get("n_total", 0)))
    hit_rows = "".join(
        f"<tr><td class='mono'>{h['candidate']}</td><td>{h['name']}</td>"
        f"<td>{h['protein']}</td><td>{h['match_type']}</td><td>{h['label']}</td>"
        f"<td>{h['host_serum']}</td>"
        f"<td><a href='https://doi.org/{h['doi']}'>{h['pmid']}</a></td></tr>"
        for h in (cv.get("hits") or []))
    miss = "".join(f"<li>{v}</li>" for v in (cv.get("missing_sequences") or {}).values())
    ttl = ("Küratörlü doğrulama seti (organizmaya özgü)" if tr
           else "Curated validation set (organism-specific)")
    return f"""<h2>{ttl}</h2>
<div class="warn"><b>{'Neden' if tr else 'Why'}:</b> {cv.get('why', '')}</div>
<p><b>{'Kapsam' if tr else 'Scope'}:</b> {cv.get('scope_note', '')}</p>
<p>{'Set' if tr else 'Set'}: <b>{cv.get('n_total', 0)}</b> {'peptit' if tr else 'peptides'}
({cv.get('n_positive', 0)} {'pozitif' if tr else 'positive'},
{cv.get('n_negative', 0)} {'negatif' if tr else 'negative'}) ·
{'sürüm' if tr else 'version'} {cv.get('version')} ({cv.get('compiled')}) ·
{'yöntem' if tr else 'method'}: <code>{cv.get('method', '')}</code></p>
<p>{'Eşleşen aday' if tr else 'Matched candidates'}:
<b>{cv.get('n_matched', 0)}</b> / {cv.get('n_candidates', 0)}</p>
<p><b>{'Payda konağa göre' if tr else 'Denominator by host'}</b> —
{'insan serumuyla doğrulanmış bir epitobu fare taramasında kaçırmak başarısızlık değildir'
 if tr else 'missing a human-serum-validated epitope in a mouse scan is not a failure'}:</p>
<table><tr><th>{'Konak' if tr else 'Host'}</th><th>{'Toplam' if tr else 'Total'}</th>
<th>{'Pozitif' if tr else 'Positive'}</th><th>{'Negatif' if tr else 'Negative'}</th></tr>
{host_rows}</table>
{"<p><b>" + ('Eşleşmeler' if tr else 'Matches') + "</b></p><table><tr><th>"
 + ('Aday' if tr else 'Candidate') + "</th><th>" + ('Epitop' if tr else 'Epitope')
 + "</th><th>" + ('Protein' if tr else 'Protein') + "</th><th>"
 + ('Tip' if tr else 'Type') + "</th><th>" + ('Etiket' if tr else 'Label')
 + "</th><th>" + ('Konak' if tr else 'Host') + "</th><th>PMID</th></tr>"
 + hit_rows + "</table>" if hit_rows else
 "<p><i>" + ('Hiçbir aday bu setle eşleşmedi. Set küçük olduğu için bu bir '
             'başarısızlık göstergesi DEĞİLDİR.' if tr else
             'No candidate matched this set. Given its small size this is NOT '
             'an indicator of failure.') + "</i></p>"}
{"<p><b>" + ('Dizisi yayınlanmadığı için sete alınamayanlar' if tr else
             'Excluded because sequences were not published') + "</b></p><ul>"
 + miss + "</ul>" if miss else ""}
<p class="mono">PMID: {', '.join(cv.get('pmids', []))}</p>"""


def _protein_label_set_section(meta: dict) -> str:
    """PROTEİN düzeyi etiket seti — Wunder 2021 bağışıklık korelatları.

    Peptit düzeyi küratörlü setin (yukarıda) protein düzeyi karşılığı. Payda
    KONAĞA göre ayrı verilir; etiketin koruma OLMADIĞI açıkça yazılır.
    """
    ps = meta.get("protein_label_set") or {}
    if not ps:
        return ""
    tr = meta.get("lang", "tr") != "en"
    ttl = ("Protein düzeyi etiket seti (bağışıklık korelatları)" if tr
           else "Protein-level label set (correlates of immunity)")
    head = (f"<h2>🧷 {ttl}</h2>"
            f"<div class='warn'><b>{'Etiketin okuması' if tr else 'Readout'}:</b> "
            f"{ps.get('readout', '')}</div>"
            f"<p class='mono'>{ps.get('citation', '')}</p>"
            f"<p>{'Set' if tr else 'Set'}: <b>{ps.get('n_features', '—')}</b> "
            f"{'array özelliği' if tr else 'array features'} / "
            f"<b>{ps.get('n_proteins', '—')}</b> "
            f"{'benzersiz protein' if tr else 'unique proteins'} · "
            f"{'koşuda bulunan locus_tag' if tr else 'locus tags in run'}: "
            f"{ps.get('n_run_locus_tags', 0)} · "
            f"{'yöntem' if tr else 'method'}: <code>{ps.get('method', '')}</code></p>")
    if ps.get("note"):
        return head + f"<p class='warn'><i>{ps['note']}</i></p>"
    rows = ""
    for h, d in (ps.get("hosts") or {}).items():
        extra = ("" if d.get("in_input") is None else
                 f"<td>{d.get('in_input')}</td><td>{d.get('kept_of_in_input')}</td>"
                 f"<td>{d.get('kept_percent_of_in_input')}</td>")
        rows += (f"<tr><td>{h}</td><td>{d.get('denominator')}</td>"
                 f"<td>{d.get('recovered')}</td><td>{d.get('recall_percent')}</td>"
                 f"<td>{d.get('unmeasured_in_host')}</td>{extra}</tr>")
    hdr = ((("Konak", "Payda", "Yakalanan", "Recall %", "O konakta ölçülmemiş")
            if tr else ("Host", "Denominator", "Recovered", "Recall %",
                        "Not measured in host")))
    hdr2 = ((("Girdide", "Elemeden geçen", "%") if tr else
             ("In input", "Kept by filters", "%"))
            if any(d.get("in_input") is not None
                   for d in (ps.get("hosts") or {}).values()) else ())
    ths = "".join(f"<th>{c}</th>" for c in (*hdr, *hdr2))
    note = ("Payda konağa göre AYRIDIR: hamsterda anlamlı bir hedefi fare "
            "taramasında kaçırmak başarısızlık değildir. 'Girdide' sütunu, "
            "'girdide hiç yoktu' durumunu 'elemede düştü' durumundan ayırır."
            if tr else
            "The denominator is host-specific: missing a hamster-significant "
            "target in a mouse scan is not a failure. The 'in input' column "
            "separates 'was never in the input' from 'dropped by the filters'.")
    return (head + f"<p>{note}</p><table><tr>{ths}</tr>{rows}</table>"
            + (f"<p><i>{ps.get('specificity_note', '')}</i></p>"
               if ps.get("specificity_note") else ""))


def _html(df: pd.DataFrame, meta: dict, peptides=None) -> str:
    lang = meta.get("lang", "tr")
    thr_rows = "".join(
        f"<tr><td>{r['step']}</td><td>{r['tool']}</td><td>{r['param']}</td>"
        f"<td>{r['value']} {r['unit']}</td><td>{'🔒' if r['hard_filter'] else '◦'}</td></tr>"
        for r in meta.get("thresholds", [])
    )
    steps_rows = "".join(
        f"<tr><td>{s['#']}</td><td>{s['adım']}</td><td>{s['durum']}</td><td>{s['not']}</td></tr>"
        for s in meta.get("plan", [])
    )
    # Atıf satırı aracın bu koşuda hazır OLMADIĞINI biliyorsa açıkça işaretlenir
    # (used is False). used None = eşlenmemiş/saf-Python yöntem → iddia yok.
    _lang_ref = "en" if meta.get("lang", "tr") == "en" else "tr"
    _absent = t(_lang_ref, "tool_status_absent")
    ref_items = "".join(
        f'<li><b>{r["tool"]}</b> ({r["step"]})'
        + ('' if r.get("used") is not False
           else f' <span style="color:#a8601a"><b>[{_absent}]</b></span>')
        + f': {r["citation"]} '
        f'<a href="{r["doi"] if r["doi"].startswith("http") else "https://doi.org/"+r["doi"]}">'
        f'{r["doi"]}</a></li>'
        for r in meta.get("citations", [])
    )
    correlate_html = _correlate_section(meta)
    curated_html = _curated_validation_section(meta)
    prot_label_html = _protein_label_set_section(meta)
    iedb_html = _iedb_section(meta, peptides)
    mev_html = _mev_section(meta)
    immuno_html = _immuno_section(meta)
    structure_html = _structure_section(meta)
    glyco_html = _glyco_section(meta)
    conservation_html = _conservation_section(meta)
    mimicry_html = _mimicry_section(meta)
    cloning_html = _cloning_section(meta)
    diagnostics_html = _diagnostics_section(meta)
    tool_status_html = _tool_status_section(meta)
    methods_html = _methods_section(meta)
    popcov = meta.get("population_coverage") or {}
    pop_html = ""
    if popcov:
        areas = popcov.get("areas", [])
        blocks = []
        for hname, he in popcov.get("hosts", {}).items():
            if he.get("note"):
                blocks.append(f'<p><b>{he.get("label", hname)}:</b> '
                              f'<i>{he["note"]}</i></p>')
                continue
            trows = ""
            for cls, lbl in (("mhc_i", "MHC-I"), ("mhc_ii", "MHC-II")):
                cov = he.get(cls)
                if not cov:
                    continue
                cells = "".join(f"<td>{cov.get(a, {}).get('coverage', '—')}%</td>"
                                if a in cov else "<td>—</td>" for a in areas)
                trows += f"<tr><td>{lbl}</td>{cells}</tr>"
            if trows:
                head = "".join(f"<th>{a}</th>" for a in areas)
                blocks.append(f'<p><b>{he.get("label", hname)}</b></p>'
                              f'<table><tr><th>{t(lang,"col_type")}</th>{head}</tr>{trows}</table>')
        note = ("" if popcov.get("available")
                else f'<p><i>{popcov.get("note", "")}</i></p>')
        pop_html = (f"<h2>{t(lang,'pop_title')}</h2>"
                    f"<p style='font-size:.85rem'>{t(lang,'pop_text')}</p>" + note + "".join(blocks))
    return f"""<!doctype html><html lang="{lang}"><head><meta charset="utf-8">
<title>{t(lang,'rep_title')}</title><style>
body{{font-family:system-ui,Arial;margin:2rem;color:#1a1a1a;max-width:1100px}}
h1{{color:#0b6}}h2{{border-bottom:2px solid #eee;padding-bottom:.3rem;margin-top:2rem}}
table{{border-collapse:collapse;width:100%;font-size:.9rem;margin:.5rem 0}}
th,td{{border:1px solid #ddd;padding:.35rem .5rem;text-align:left}}
th{{background:#f4f4f4}}code{{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px}}
.warn{{background:#fff4e5;border-left:4px solid #b26a00;padding:.6rem;margin:1rem 0}}
.mono{{font-family:monospace;font-size:.72rem;color:#555;word-break:break-all}}
td.pass{{color:#0a7;font-weight:600}}td.fail{{color:#c33;font-weight:600}}td.na{{color:#999}}
span.pass{{color:#0a7;font-weight:700}}span.fail{{color:#c33;font-weight:600}}
span.warn{{color:#b26a00;font-weight:600}}
tr.sel{{background:#e7f7ee}}tr.sel td{{border-color:#b7e6ca}}
details{{margin:.4rem 0;border:1px solid #e5e5e5;border-radius:5px;padding:.3rem .6rem}}
summary{{cursor:pointer;font-size:.95rem}}</style></head><body>
<h1>🧬 {t(lang,'rep_title')}</h1>
<p>{t(lang,'rep_generated')}: {meta.get('timestamp','')} · {t(lang,'rep_input')}: <code>{meta.get('input','')}</code> ·
{t(lang,'rep_profile')}: <code>{meta.get('profile','')}</code></p>
<div class="warn"><b>Not/Note:</b> {t(lang,'rep_disclaimer')}</div>

{correlate_html}

<h2>{t(lang,'rep_summary')}</h2>
<table><tr><th>{t(lang,'rep_sum_input')}</th><th>{t(lang,'rep_sum_discovery')}</th><th>{t(lang,'rep_sum_funnel')}</th>
<th>{t(lang,'rep_sum_epitope')}</th><th>{t(lang,'rep_sum_survivors')}</th></tr>
<tr><td>{meta.get('n_input','?')}</td><td>{meta.get('n_discovery','?')}</td>
<td>{meta.get('n_funnel','?')}</td><td>{meta.get('n_epitope','?')}</td>
<td>{meta.get('n_survivors','?')}</td></tr></table>

{methods_html}

<h2>{t(lang,'rep_top15')}</h2>
{df.head(15).to_html(index=False)}

<h2>{t(lang,'rep_bytype')}</h2>
{_type_tables_html(peptides or [], meta)}

{iedb_html}
{curated_html}
{prot_label_html}

{mev_html}
{structure_html}
{cloning_html}
{tool_status_html}
{diagnostics_html}
{conservation_html}
{glyco_html}
{mimicry_html}
{immuno_html}

{pop_html}

<h2>{t(lang,'rep_thresholds')}</h2>
<table><tr><th>{t(lang,'col_step')}</th><th>{t(lang,'col_tool')}</th><th>{t(lang,'col_param')}</th><th>{t(lang,'col_value')}</th><th>{t(lang,'col_type')}</th></tr>{thr_rows}</table>

<h2>{t(lang,'rep_plan')}</h2>
<table><tr><th>#</th><th>{t(lang,'col_step')}</th><th>{t(lang,'col_status')}</th><th>{t(lang,'col_note')}</th></tr>{steps_rows}</table>

<!-- REFERANSLAR EN SONDA (Ali kuralı 2026-11-22) — PDF ile aynı sıra. -->
<h2>Kaynaklar / atıflar</h2>
<ol style="font-size:.8rem">{ref_items}</ol>
</body></html>"""
