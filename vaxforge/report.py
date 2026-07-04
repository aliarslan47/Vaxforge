"""Adım 8 — Rapor + makine-okur veri paketi + HTML panosu.

Üretilenler (outputs/<run>/):
  candidates.csv     — sıralı aday peptitler + metrikler
  top_peptides.fasta — en iyi adaylar
  construct.gb       — mRNA konstrüktü (GenBank)
  run.json           — tüm parametreler/eşikler/özetler (tekrarlanabilirlik)
  report.html        — insan-okur panosu
PDF sonraki sürümde (reportlab/weasyprint).
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .mrna import Construct
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
            "allergen": p.metrics.get("allergen"),
            "toxicity": p.metrics.get("toxicity"),
            "method": p.methods.get("mhc_score"),
            "passed": p.passed,
        })
    return pd.DataFrame(rows)


def write_package(outdir: Path, peptides: list[Peptide], construct: Construct,
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

    # GenBank mRNA
    rec = SeqRecord(Seq(construct.mrna), id="VaxForge_mRNA",
                    description="çok-epitoplu mRNA aşı konstrüktü (in silico)")
    rec.annotations["molecule_type"] = "mRNA"
    rec.annotations["topology"] = "linear"
    p_gb = outdir / "construct.gb"
    SeqIO.write(rec, p_gb, "genbank")
    paths["genbank"] = p_gb

    # JSON (tam koşum)
    run_meta = dict(run_meta)
    run_meta["construct"] = {"protein": construct.protein, "mrna": construct.mrna,
                             "parts": construct.parts, "metrics": construct.metrics,
                             "methods": construct.methods}
    run_meta["candidates"] = df.to_dict(orient="records")
    p_json = outdir / "run.json"
    p_json.write_text(json.dumps(run_meta, indent=2, ensure_ascii=False, default=str))
    paths["json"] = p_json

    # HTML panosu
    p_html = outdir / "report.html"
    p_html.write_text(_html(df, construct, run_meta), encoding="utf-8")
    paths["html"] = p_html

    # PDF (yayın-tarzı) — reportlab varsa
    try:
        from . import report_pdf
        paths["pdf"] = report_pdf.build(outdir, peptides, construct, run_meta)
    except Exception as e:
        (outdir / "pdf_error.txt").write_text(f"PDF üretilemedi: {e}")
    return paths


def _html(df: pd.DataFrame, construct: Construct, meta: dict) -> str:
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
    cm = construct.metrics
    return f"""<!doctype html><html lang="tr"><head><meta charset="utf-8">
<title>VaxForge Raporu</title><style>
body{{font-family:system-ui,Arial;margin:2rem;color:#1a1a1a;max-width:1100px}}
h1{{color:#0b6}}h2{{border-bottom:2px solid #eee;padding-bottom:.3rem;margin-top:2rem}}
table{{border-collapse:collapse;width:100%;font-size:.9rem;margin:.5rem 0}}
th,td{{border:1px solid #ddd;padding:.35rem .5rem;text-align:left}}
th{{background:#f4f4f4}}code{{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px}}
.warn{{background:#fff4e5;border-left:4px solid #b26a00;padding:.6rem;margin:1rem 0}}
.mono{{font-family:monospace;font-size:.75rem;word-break:break-all}}</style></head><body>
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

<h2>En iyi 15 aday peptit</h2>
{df.head(15).to_html(index=False)}

<h2>mRNA konstrüktü</h2>
<p>Yapı: {' + '.join(construct.parts)}</p>
<table><tr><th>Protein uz.</th><th>mRNA uz.</th><th>GC%</th><th>CAI (insan)</th>
<th>Kararsızlık</th><th>pI</th><th>CTL/HTL/B</th></tr>
<tr><td>{cm.get('protein_len')}</td><td>{cm.get('mrna_len')}</td><td>{cm.get('gc_percent')}</td>
<td>{cm.get('cai_human')}</td><td>{cm.get('instability')}</td><td>{cm.get('pI')}</td>
<td>{cm.get('n_ctl')}/{cm.get('n_htl')}/{cm.get('n_bcell')}</td></tr></table>
<p class="mono"><b>mRNA:</b> {construct.mrna}</p>

<h2>Kullanılan eşikler (tekrarlanabilirlik)</h2>
<table><tr><th>Adım</th><th>Araç</th><th>Parametre</th><th>Değer</th><th>Tip</th></tr>{thr_rows}</table>

<h2>Çalıştırılan plan</h2>
<table><tr><th>#</th><th>Adım</th><th>Durum</th><th>Not</th></tr>{steps_rows}</table>

<h2>Kaynaklar / atıflar</h2>
<ol style="font-size:.8rem">{ref_items}</ol>
</body></html>"""
