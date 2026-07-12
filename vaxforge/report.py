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

from .i18n import t
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
    titles = {
        "MHC-I": f"🟦 {t(lang,'tt_ctl')} <span class='mono'>({t(lang,'tt_rank_asc')})</span>",
        "MHC-II": f"🟨 {t(lang,'tt_htl')} <span class='mono'>({t(lang,'tt_rank_asc')})</span>",
        "B": f"🟩 {t(lang,'tt_bcell')} <span class='mono'>({t(lang,'tt_bepi_desc')})</span>"}

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
        trs = f"<tr><td colspan='8'><i>{t(lang,'iedb_none')}</i></td></tr>"
    bm = im.get("benchmark") or {}
    bm_html = ""
    if bm and bm.get("recall") is not None:
        k = bm.get("k", 8)
        bm_html = (
            f"<h3>{t(lang,'iedb_val_title')}</h3>"
            f"<p style='font-size:.85rem'>{t(lang,'iedb_val_text').format(k=k)}</p>"
            f"<table><tr><th>{t(lang,'col_known')}</th><th>{t(lang,'col_captured')}</th>"
            f"<th>{t(lang,'col_recall')}</th><th>{t(lang,'col_pred')}</th><th>{t(lang,'col_pred_match')}</th>"
            f"<th>{t(lang,'col_matchrate')}</th></tr>"
            f"<tr><td>{bm['n_known']}</td><td>{bm['n_known_hit']}</td>"
            f"<td><b>{round(bm['recall']*100,1)}%</b></td><td>{bm['n_pred']}</td>"
            f"<td>{bm['n_pred_matched']}</td>"
            f"<td><b>{round((bm['precision_like'] or 0)*100,1)}%</b></td></tr></table>"
        )
    return (
        f"<h2>{t(lang,'iedb_title')}</h2>"
        f"<p style='font-size:.85rem'>{t(lang,'iedb_intro').format(src=im.get('source',''), n=im.get('n_matched',0), tot=im.get('n_candidates','?'))}</p>"
        f"<table><tr><th>{t(lang,'col_order')}</th><th>{t(lang,'col_candidate')}</th><th>{t(lang,'col_type')}</th><th>{t(lang,'col_match')}</th>"
        f"<th>{t(lang,'col_iedb_epi')}</th><th>{t(lang,'col_organism')}</th><th>{t(lang,'col_antigen')}</th><th>{t(lang,'col_reference')}</th></tr>"
        f"{trs}</table>{bm_html}"
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

<h2>{t(lang,'rep_top15')}</h2>
{df.head(15).to_html(index=False)}

<h2>{t(lang,'rep_bytype')}</h2>
{_type_tables_html(peptides or [], meta)}

<h2>{t(lang,'rep_thresholds')}</h2>
<table><tr><th>{t(lang,'col_step')}</th><th>{t(lang,'col_tool')}</th><th>{t(lang,'col_param')}</th><th>{t(lang,'col_value')}</th><th>{t(lang,'col_type')}</th></tr>{thr_rows}</table>

<h2>{t(lang,'rep_plan')}</h2>
<table><tr><th>#</th><th>{t(lang,'col_step')}</th><th>{t(lang,'col_status')}</th><th>{t(lang,'col_note')}</th></tr>{steps_rows}</table>

{iedb_html}

{pop_html}

<h2>Kaynaklar / atıflar</h2>
<ol style="font-size:.8rem">{ref_items}</ol>
</body></html>"""
