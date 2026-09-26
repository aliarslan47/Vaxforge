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
from datetime import datetime
from pathlib import Path

import pandas as pd

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
            "anchor_residues": ";".join(f"{k}={v}" for k, v in
                                        (p.metrics.get("anchor_residues") or {}).items()),
            "allele_anchor_motif": ";".join(f"{k}∈{{{','.join(v)}}}" for k, v in
                                            (p.metrics.get("allele_anchor_motif") or {}).items()),
            "anchor_match": p.metrics.get("anchor_match"),
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
            "method": method_label(p.methods.get("mhc_score"), lang),
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
    p_fa = outdir / "top_peptides.fasta"
    with p_fa.open("w") as fh:
        for i, p in enumerate(peptides[:20], 1):
            g = p.metrics.get("gene") or ""
            lt = p.metrics.get("locus_tag") or ""
            fh.write(f">cand{i}|{p.kind}|score={p.candidacy}|cds={p.parent}"
                     f"{('|gene='+g) if g else ''}{('|locus='+lt) if lt else ''}\n{p.seq}\n")
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

    # JSON (tam koşum)
    run_meta = dict(run_meta)
    run_meta["candidates"] = df.to_dict(orient="records")
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
        return '<span class="pass">✔</span>' if b else '<span class="fail">✗</span>'

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
            base = (f'<td>{star}</td><td class="mono"><b>{r["epitope"]}</b></td>'
                    f'<td class="mono">{r["source"]}</td><td>{r["pos"]}</td>')
            if kind == "B":
                trs += (f'<tr{cls}>{base}<td>{r["length"]}</td><td>{nz(r["bcell"])}</td>'
                        f'<td>{nz(r["antigenicity"])}</td><td>{ok(r["allergen_ok"])}</td>'
                        f'<td>{ok(r["toxic_ok"])}</td><td>{lit}</td></tr>')
            elif kind == "MHC-I":
                trs += (f'<tr{cls}>{base}<td>{nz(r["rank"])}</td><td class="mono">{r["allele"]}</td>'
                        f'<td>{nz(r["antigenicity"])}</td><td>{ok(r["allergen_ok"])}</td>'
                        f'<td>{ok(r["toxic_ok"])}</td><td>{nz(r["immunogenicity"])}</td>'
                        f'<td>{nz(r["processing"])}</td><td>–</td><td>{lit}</td></tr>')
            else:
                trs += (f'<tr{cls}>{base}<td>{nz(r["rank"])}</td><td class="mono">{r["allele"]}</td>'
                        f'<td>{nz(r["antigenicity"])}</td><td>{ok(r["allergen_ok"])}</td>'
                        f'<td>{ok(r["toxic_ok"])}</td><td>{ok(r["ifn_ok"])}</td><td>–</td><td>{lit}</td></tr>')
        out.append(f"<table>{head}{trs}</table>")
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
        ("Allergenicity", yn(al["allergenic"], "Allergen", "Non-allergen") + f" · {al['method']}"),
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
    return (
        f"<h2>🧩 {ttl}</h2>"
        f"<p style='font-size:.85rem'>{comp} · {len_lbl}: <b>{len(mv['sequence'])} aa</b> "
        f"(adjuvan: {mv['adjuvant']})</p>"
        f"<details><summary>{seq_lbl}</summary>"
        f"<p class='mono'>{mv['sequence']}</p></details>"
        f"<table><tr><th>{t(lang,'col_param') if False else ('Property' if is_en else 'Özellik')}</th>"
        f"<th>{'Value' if is_en else 'Değer'}</th></tr>{trows}</table>"
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
    hdr = (["Host", "Verdict", "Score", "MHC-I", "MHC-II", "Clearance", "Profile", "Ig"] if is_en
           else ["Konak", "Karar", "Skor", "MHC-I", "MHC-II", "Temizlik", "Profil", "Ig"])
    trs = ""
    for ph in imm["per_host"]:
        pr = ph.get("profile") or {}
        prof = f"{pr.get('polarization', '—')} ({round(pr.get('th1', 0) * 100)}/{round(pr.get('th2', 0) * 100)})"
        trs += ("<tr>"
                f"<td>{ph.get('host_label', ph.get('host', ''))}</td>"
                f"<td>{vmap.get(ph.get('verdict', ''), ph.get('verdict', ''))}</td>"
                f"<td>{round(ph.get('score', 0))}/100</td>"
                f"<td>{ph.get('n_mhci', '')}</td><td>{ph.get('n_mhcii', '')}</td>"
                f"<td>{round(ph.get('cleared_fraction', 0) * 100)}%</td>"
                f"<td>{prof}</td><td>{pr.get('dominant_isotype') or '—'}</td></tr>")
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
    rows = [
        (("Source" if is_en else "Kaynak"), st.get("source", "—")),
        (("Residues" if is_en else "Kalıntı"), st.get("n_res", "—")),
        ("pLDDT (mean/min)", f"{st.get('plddt_mean','—')} / {st.get('plddt_min','—')}"),
        (("Surface %" if is_en else "Yüzey %"), f"{surf.get('percent_surface','—')}%"),
        (("Conformational B patches" if is_en else "Konformasyonel B yaması"),
         cb.get("n_patches", "—")),
    ]
    trs = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    patch_html = ""
    for p in (cb.get("patches") or [])[:6]:
        patch_html += (f"<tr><td>{p['start']}–{p['end']}</td><td>{p['length']}</td>"
                       f"<td class='mono'>{p['peptide']}</td><td>{p['mean_score']}</td></tr>")
    if patch_html:
        ph = ("<p style='margin-top:.6rem'><b>Surface B-cell patches</b></p>" if is_en
              else "<p style='margin-top:.6rem'><b>Yüzey B-hücre yamaları</b></p>")
        hh = (["Range", "Len", "Peptide", "Score"] if is_en else ["Aralık", "Uzunluk", "Peptit", "Skor"])
        patch_html = ph + "<table><tr>" + "".join(f"<th>{h}</th>" for h in hh) + f"</tr>{patch_html}</table>"
    return (
        f"<h2>🧬 {ttl}</h2>"
        f"<p style='font-size:.8rem;background:#eef6ff;border:1px solid #9cc3e6;padding:.5rem;border-radius:6px'>ℹ {disc}</p>"
        f"<table>{trs}</table>{patch_html}"
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
    rows = [
        (("Host" if is_en else "Konak"), cl.get("host", "—")),
        (("Protein length (aa)" if is_en else "Protein (aa)"), cl.get("protein_len", "—")),
        (("DNA length (bp)" if is_en else "DNA (bp)"), cl.get("dna_len_bp", "—")),
        ("CAI", cl.get("cai", "—")),
        ("GC %", f"{cl.get('gc_percent', '—')}%"),
        (("Internal restriction sites" if is_en else "İç restriksiyon bölgeleri"), sites_str),
        (("Suggested vector" if is_en else "Önerilen vektör"), cl.get("suggested_vector", "—")),
        (("Tag / strategy" if is_en else "Etiket / strateji"), cl.get("tag", "—")),
    ]
    trs = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    dna = cl.get("dna", "")
    dna_prev = (dna[:120] + "…") if len(dna) > 120 else dna
    return (
        f"<h2>🧫 {ttl}</h2>"
        f"<p style='font-size:.8rem;background:#eef6ff;border:1px solid #9cc3e6;padding:.5rem;border-radius:6px'>ℹ {disc}</p>"
        f"<table>{trs}</table>"
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
    return (
        f"<h2>🧪 {t(lang, 'rep_methods')}</h2>"
        f"<p style='font-size:.85rem'>{t(lang, 'rep_methods_sub')}</p>"
        f"<table><tr><th>{t(lang, 'meth_step')}</th><th>{t(lang, 'meth_method')}</th></tr>{trs}</table>"
    )


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
    ref_items = "".join(
        f'<li><b>{r["tool"]}</b> ({r["step"]}): {r["citation"]} '
        f'<a href="{r["doi"] if r["doi"].startswith("http") else "https://doi.org/"+r["doi"]}">'
        f'{r["doi"]}</a></li>'
        for r in meta.get("citations", [])
    )
    iedb_html = _iedb_section(meta, peptides)
    mev_html = _mev_section(meta)
    immuno_html = _immuno_section(meta)
    structure_html = _structure_section(meta)
    glyco_html = _glyco_section(meta)
    mimicry_html = _mimicry_section(meta)
    cloning_html = _cloning_section(meta)
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
span.pass{{color:#0a7;font-weight:700}}span.fail{{color:#c33}}
tr.sel{{background:#e7f7ee}}tr.sel td{{border-color:#b7e6ca}}
details{{margin:.4rem 0;border:1px solid #e5e5e5;border-radius:5px;padding:.3rem .6rem}}
summary{{cursor:pointer;font-size:.95rem}}</style></head><body>
<h1>🧬 {t(lang,'rep_title')}</h1>
<p>{t(lang,'rep_generated')}: {meta.get('timestamp','')} · {t(lang,'rep_input')}: <code>{meta.get('input','')}</code> ·
{t(lang,'rep_profile')}: <code>{meta.get('profile','')}</code></p>
<div class="warn"><b>Not/Note:</b> {t(lang,'rep_disclaimer')}</div>

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

{mev_html}
{structure_html}
{cloning_html}
{glyco_html}
{mimicry_html}
{immuno_html}

{iedb_html}

{pop_html}

<h2>Kaynaklar / atıflar</h2>
<ol style="font-size:.8rem">{ref_items}</ol>

<h2>{t(lang,'rep_thresholds')}</h2>
<table><tr><th>{t(lang,'col_step')}</th><th>{t(lang,'col_tool')}</th><th>{t(lang,'col_param')}</th><th>{t(lang,'col_value')}</th><th>{t(lang,'col_type')}</th></tr>{thr_rows}</table>

<h2>{t(lang,'rep_plan')}</h2>
<table><tr><th>#</th><th>{t(lang,'col_step')}</th><th>{t(lang,'col_status')}</th><th>{t(lang,'col_note')}</th></tr>{steps_rows}</table>
</body></html>"""
