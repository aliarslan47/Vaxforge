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

from .models import Peptide


def _candidates_df(peptides: list[Peptide]) -> pd.DataFrame:
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
            "anchor_residues": ";".join(f"{k}={v}" for k, v in
                                        (p.metrics.get("anchor_residues") or {}).items()),
            "allele_anchor_motif": ";".join(f"{k}∈{{{','.join(v)}}}" for k, v in
                                            (p.metrics.get("allele_anchor_motif") or {}).items()),
            "anchor_match": p.metrics.get("anchor_match"),
            "allergen": p.metrics.get("allergen"),
            "toxicity": p.metrics.get("toxicity"),
            "iedb_match": _iedb_cell(p, "type"),
            "iedb_epitope": _iedb_cell(p, "epitope"),
            "iedb_organism": _iedb_cell(p, "organism"),
            "iedb_pmids": _iedb_cell(p, "pmids"),
            "method": p.methods.get("mhc_score"),
            "passed": p.passed,
        })
    return pd.DataFrame(rows)


def _iedb_cell(p: Peptide, field: str):
    ie = p.metrics.get("iedb") or {}
    if ie.get("matched") is not True:
        return "" if field != "type" else ("yok" if ie.get("matched") is False else "n/a")
    if field == "type":
        return ie.get("match_type", "eşleşme")
    if field == "epitope":
        return ie.get("epitope_seq", "")
    if field == "organism":
        return (ie.get("organisms") or [""])[0]
    if field == "pmids":
        return ";".join(ie.get("pmids") or [])
    return ""


def write_package(outdir: Path, peptides: list[Peptide],
                  run_meta: dict) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {}

    # CSV
    df = _candidates_df(peptides)
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


def _candidate_blocks(peptides, meta) -> str:
    """Öne çıkan adaylar (en iyi 15 + literatürde eşleşen) için araç sonuçları +
    GEÇTİ/GEÇEMEDİ tablosu. Tam liste candidates_full.xlsx'te."""
    from . import evaluate
    thr = evaluate.thr_lookup(meta)
    subset = evaluate.report_subset(peptides, top_n=15)
    intro = (f'<p style="font-size:.85rem">En iyi 15 aday + literatürde (IEDB) eşleşen '
             f'adaylar gösterilir ({len(subset)} / {len(peptides)} aday). Tüm adayların '
             f'tam araç-sonucu tablosu <code>candidates_full.xlsx</code> dosyasındadır.</p>')
    blocks = [intro]
    for n, (i, p, reason) in enumerate(subset, 1):
        rows = evaluate.candidate_rows(p, thr)
        trs = ""
        for r in rows:
            cls = ("pass" if r["status"] == evaluate.PASS else
                   "fail" if r["status"] == evaluate.FAIL else "na")
            hard = " 🔒" if r["hard"] else ""
            trs += (f'<tr><td>{r["tool"]}{hard}</td><td>{r["value"]}</td>'
                    f'<td>{r["cutoff"]}</td><td class="{cls}">{r["status"]}</td>'
                    f'<td class="mono">{r["method"]}</td></tr>')
        gene = p.metrics.get("gene") or "—"
        tag = ' · 📖 <b>literatürde</b>' if reason == "literatür" else ""
        blocks.append(
            f'<details {"open" if n <= 3 else ""}><summary><b>#{i} {p.seq}</b> '
            f'— {p.kind} — adaylık <b>{p.candidacy}</b> '
            f'(kaynak: {p.parent}, gen: {gene}){tag}</summary>'
            f'<table><tr><th>Araç</th><th>Sonuç</th><th>Eşik (referans)</th>'
            f'<th>Durum</th><th>Yöntem</th></tr>{trs}</table></details>')
    return "".join(blocks)


def _iedb_section(meta: dict, peptides) -> str:
    """IEDB literatür/bilinen-epitop eşleşmesi + (varsa) recall benchmark bölümü."""
    im = meta.get("iedb_match")
    if not im:
        return ""
    if not im.get("available"):
        return ('<h2>IEDB literatür/bilinen-epitop taraması</h2>'
                f'<div class="warn">{im.get("note", "IEDB taraması yapılamadı")}</div>')
    peptides = peptides or []
    # eşleşen adaylar tablosu
    trs = ""
    for i, p in enumerate(peptides, 1):
        ie = p.metrics.get("iedb") or {}
        if ie.get("matched") is not True:
            continue
        org = "; ".join((ie.get("organisms") or [])[:2]) or "—"
        ag = "; ".join((ie.get("antigens") or [])[:2]) or "—"
        pmids = ie.get("pmids") or []
        plinks = ", ".join(
            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pm}/">{pm}</a>' for pm in pmids[:5]
        ) or (f'<a href="https://www.iedb.org/epitope/{ie.get("iedb_id")}">IEDB {ie.get("iedb_id")}</a>'
              if ie.get("iedb_id") else "—")
        trs += (f"<tr><td>#{i}</td><td class='mono'>{p.seq}</td><td>{p.kind}</td>"
                f"<td>{ie.get('match_type','')}</td><td class='mono'>{ie.get('epitope_seq','')}</td>"
                f"<td>{org}</td><td>{ag}</td><td>{plinks}</td></tr>")
    if not trs:
        trs = "<tr><td colspan='8'><i>Hiçbir aday bilinen IEDB epitobuyla eşleşmedi.</i></td></tr>"
    bm = im.get("benchmark") or {}
    bm_html = ""
    if bm and bm.get("recall") is not None:
        k = bm.get("k", 8)
        bm_html = (
            "<h3>Validasyon — bilinen epitop recall'ü</h3>"
            "<p style='font-size:.85rem'>Bu organizma için IEDB'de deneysel doğrulanmış "
            "lineer epitoplar 'ground truth' alınır; pipeline'ın tahminleri bunlarla "
            f"örtüşme (≥{k} aa ortak çekirdek / içerme / exact) üzerinden değerlendirilir. "
            "NOT: pipeline seçici olarak KISA bir öncelik listesi üretir (protein başına "
            "sınırlı peptit); bu yüzden tüm epitop kataloğuna karşı recall doğası gereği "
            "düşüktür — asıl anlamlı ölçüt, adayların ne kadarının deneysel doğrulanmış "
            "olduğudur (eşleşme oranı).</p>"
            "<table><tr><th>Bilinen epitop (benzersiz)</th><th>Yakalanan</th>"
            "<th>Recall</th><th>Tahmin (aday)</th><th>Bilinene eşleşen aday</th>"
            "<th>Eşleşme oranı</th></tr>"
            f"<tr><td>{bm['n_known']}</td><td>{bm['n_known_hit']}</td>"
            f"<td><b>{round(bm['recall']*100,1)}%</b></td><td>{bm['n_pred']}</td>"
            f"<td>{bm['n_pred_matched']}</td>"
            f"<td><b>{round((bm['precision_like'] or 0)*100,1)}%</b></td></tr></table>"
        )
    return (
        "<h2>IEDB literatür/bilinen-epitop taraması</h2>"
        f"<p style='font-size:.85rem'>Kaynak: <code>{im.get('source','')}</code> · "
        f"eşleşen aday: <b>{im.get('n_matched',0)}/{im.get('n_candidates','?')}</b>. "
        "Eşleşme, adayın deneysel doğrulanmış bir epitopla (exact/içerme/ortak çekirdek) "
        "örtüştüğünü gösterir — güçlü pozitif kontrol sinyali. Bu adım adaylık puanını "
        "<b>değiştirmez</b> (salt yorumlama).</p>"
        "<table><tr><th>Sıra</th><th>Aday peptit</th><th>Tip</th><th>Eşleşme</th>"
        "<th>IEDB epitobu</th><th>Organizma</th><th>Antijen</th><th>Referans</th></tr>"
        f"{trs}</table>{bm_html}"
    )


def _html(df: pd.DataFrame, meta: dict, peptides=None) -> str:
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
                              f'<table><tr><th>Sınıf</th>{head}</tr>{trows}</table>')
        note = ("" if popcov.get("available")
                else f'<p><i>{popcov.get("note", "")}</i></p>')
        pop_html = ("<h2>Popülasyon kapsamı (IEDB HLA frekansları)</h2>"
                    "<p style='font-size:.85rem'>Bir bireyin en az bir epitop-bağlayan "
                    "allele sahip olma olasılığı (%). Gerçek frekans verisi yalnız insan "
                    "HLA için mevcuttur.</p>" + note + "".join(blocks))
    return f"""<!doctype html><html lang="tr"><head><meta charset="utf-8">
<title>VaxForge Raporu</title><style>
body{{font-family:system-ui,Arial;margin:2rem;color:#1a1a1a;max-width:1100px}}
h1{{color:#0b6}}h2{{border-bottom:2px solid #eee;padding-bottom:.3rem;margin-top:2rem}}
table{{border-collapse:collapse;width:100%;font-size:.9rem;margin:.5rem 0}}
th,td{{border:1px solid #ddd;padding:.35rem .5rem;text-align:left}}
th{{background:#f4f4f4}}code{{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px}}
.warn{{background:#fff4e5;border-left:4px solid #b26a00;padding:.6rem;margin:1rem 0}}
.mono{{font-family:monospace;font-size:.72rem;color:#555;word-break:break-all}}
td.pass{{color:#0a7;font-weight:600}}td.fail{{color:#c33;font-weight:600}}td.na{{color:#999}}
details{{margin:.4rem 0;border:1px solid #e5e5e5;border-radius:5px;padding:.3rem .6rem}}
summary{{cursor:pointer;font-size:.95rem}}</style></head><body>
<h1>🧬 VaxForge — Aşı Adayı Raporu</h1>
<p>Oluşturma: {meta.get('timestamp','')} · Girdi: <code>{meta.get('input','')}</code> ·
Profil: <code>{meta.get('profile','')}</code></p>
<div class="warn"><b>Not:</b> Bu prototip, harici araçlar takılı olmadığında
saf-Python <b>heuristik/proxy</b> yöntemleri kullanır (VaxiJen, NetMHCpan, AllerTOP,
AlphaFold vb. yerine). Skorlar gerçek araç çıktısı değildir; yöntem etiketleri
JSON'da <code>methods</code> altındadır.</div>

<h2>Özet</h2>
<table><tr><th>Girdi proteini</th><th>Keşif sonrası</th><th>Huni sonrası</th>
<th>Epitop</th><th>Sağ kalan aday</th></tr>
<tr><td>{meta.get('n_input','?')}</td><td>{meta.get('n_discovery','?')}</td>
<td>{meta.get('n_funnel','?')}</td><td>{meta.get('n_epitope','?')}</td>
<td>{meta.get('n_survivors','?')}</td></tr></table>

<h2>En iyi 15 aday peptit (özet)</h2>
{df.head(15).to_html(index=False)}

<h2>Aday peptitler — TÜM araç sonuçları (en iyiden en kötüye)</h2>
<p style="font-size:.82rem">Her aday için çalışan tüm araçların çıktısı ve referans eşiğe göre
durumu. 🔒 = sert filtre (geçemeyen elenir); eşiksiz satırlar yorum amaçlıdır.</p>
{_candidate_blocks(peptides or [], meta)}

<h2>Kullanılan eşikler (tekrarlanabilirlik)</h2>
<table><tr><th>Adım</th><th>Araç</th><th>Parametre</th><th>Değer</th><th>Tip</th></tr>{thr_rows}</table>

<h2>Çalıştırılan plan</h2>
<table><tr><th>#</th><th>Adım</th><th>Durum</th><th>Not</th></tr>{steps_rows}</table>

{iedb_html}

{pop_html}

<h2>Kaynaklar / atıflar</h2>
<ol style="font-size:.8rem">{ref_items}</ol>
</body></html>"""
