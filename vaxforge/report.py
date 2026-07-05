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
            "method": p.methods.get("mhc_score"),
            "passed": p.passed,
        })
    return pd.DataFrame(rows)


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

    # JSON (tam koşum)
    run_meta = dict(run_meta)
    run_meta["candidates"] = df.to_dict(orient="records")
    p_json = outdir / "run.json"
    p_json.write_text(json.dumps(run_meta, indent=2, ensure_ascii=False, default=str))
    paths["json"] = p_json

    # HTML panosu
    p_html = outdir / "report.html"
    p_html.write_text(_html(df, run_meta, peptides), encoding="utf-8")
    paths["html"] = p_html

    # PDF (yayın-tarzı) — reportlab varsa
    try:
        from . import report_pdf
        paths["pdf"] = report_pdf.build(outdir, peptides, run_meta)
    except Exception as e:
        (outdir / "pdf_error.txt").write_text(f"PDF üretilemedi: {e}")
    return paths


def _candidate_blocks(peptides, meta) -> str:
    """Her aday için (best→worst) tüm araç sonuçları + GEÇTİ/GEÇEMEDİ tablosu."""
    from . import evaluate
    thr = evaluate.thr_lookup(meta)
    blocks = []
    for i, p in enumerate(peptides, 1):
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
        blocks.append(
            f'<details {"open" if i <= 3 else ""}><summary><b>#{i} {p.seq}</b> '
            f'— {p.kind} — adaylık <b>{p.candidacy}</b> '
            f'(kaynak: {p.parent}, gen: {gene})</summary>'
            f'<table><tr><th>Araç</th><th>Sonuç</th><th>Eşik (referans)</th>'
            f'<th>Durum</th><th>Yöntem</th></tr>{trs}</table></details>')
    return "".join(blocks)


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

{pop_html}

<h2>Kaynaklar / atıflar</h2>
<ol style="font-size:.8rem">{ref_items}</ol>
</body></html>"""
