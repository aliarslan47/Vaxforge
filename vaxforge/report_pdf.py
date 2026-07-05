"""Yayın-tarzı PDF rapor (reportlab + matplotlib).

Bilimsel çekirdek artık gerçek araçlarla çalıştığından, sunum/tez için düzgün bir
PDF üretir: özet, aday grafiği, en iyi peptitler, kullanılan eşikler,
yöntemler/araçlar ve referanslar. Türkçe için DejaVuSans fontu kaydedilir.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, Paragraph, SimpleDocTemplate, Spacer,
                                Table, TableStyle)

from . import citations

_FONT = "DejaVu"
_FONT_B = "DejaVu-Bold"


def _register_font() -> bool:
    base = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    try:
        pdfmetrics.registerFont(TTFont(_FONT, str(base / "DejaVuSans.ttf")))
        pdfmetrics.registerFont(TTFont(_FONT_B, str(base / "DejaVuSans-Bold.ttf")))
        return True
    except Exception:
        return False


def _chart(peptides, path: Path) -> bool:
    top = [p for p in peptides if p.passed][:12]
    if not top:
        return False
    colmap = {"B": "#2e7d32", "MHC-I": "#1565c0", "MHC-II": "#b26a00"}
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    labels = [f"{p.seq[:11]}" for p in top]
    vals = [p.candidacy for p in top]
    cols = [colmap.get(p.kind, "#777") for p in top]
    ax.bar(range(len(top)), vals, color=cols)
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Adaylık puanı")
    ax.set_title("En iyi aday peptitler (renk = tip)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colmap.values()]
    ax.legend(handles, colmap.keys(), fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def build(outdir: Path, peptides, meta: dict) -> Path:
    outdir = Path(outdir)
    _register_font()
    out = outdir / "report.pdf"
    styles = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=styles["Normal"], fontName=_FONT, fontSize=9, leading=12)
    h1 = ParagraphStyle("h1", parent=styles["Heading1"], fontName=_FONT_B, fontSize=16, textColor=colors.HexColor("#0b6b4f"))
    h2 = ParagraphStyle("h2", parent=styles["Heading2"], fontName=_FONT_B, fontSize=11, spaceBefore=10)
    small = ParagraphStyle("small", parent=body, fontSize=7, textColor=colors.grey)
    mono = ParagraphStyle("mono", parent=body, fontName="Courier", fontSize=6, leading=7)

    doc = SimpleDocTemplate(str(out), pagesize=A4, topMargin=1.5 * cm,
                            bottomMargin=1.5 * cm, leftMargin=1.6 * cm, rightMargin=1.6 * cm)
    el = []

    def tbl(data, widths=None, header=True):
        t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
        style = [("FONTNAME", (0, 0), (-1, -1), _FONT), ("FONTSIZE", (0, 0), (-1, -1), 7),
                 ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cccccc")),
                 ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f6f6")])]
        if header:
            style += [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0b6b4f")),
                      ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("FONTNAME", (0, 0), (-1, 0), _FONT_B)]
        t.setStyle(TableStyle(style))
        return t

    # -- Başlık
    el.append(Paragraph("VaxForge — Aşı Adayı Raporu", h1))
    hosts = ", ".join(h["label"] for h in meta.get("hosts", [])) or "—"
    el.append(Paragraph(f"Oluşturma: {meta.get('timestamp','')} &nbsp;·&nbsp; Girdi: "
                        f"<b>{meta.get('input','')}</b> &nbsp;·&nbsp; Patojen profili: "
                        f"<b>{meta.get('profile','')}</b> &nbsp;·&nbsp; Konak(lar): <b>{hosts}</b>", body))
    el.append(Spacer(1, 6))
    el.append(Paragraph("Bu rapor VaxForge in silico reverse vaccinology hattı ile üretilmiştir. "
                        "Bilimsel adımlar gerçek araçlarla koşulmuştur; kullanılan araçlar ve eşikler aşağıda listelenir.", small))

    # -- Özet
    el.append(Paragraph("1. Özet", h2))
    el.append(tbl([["Girdi proteini", "Keşif sonrası", "Huni sonrası", "Epitop", "Sağ kalan aday"],
                   [meta.get("n_input", "?"), meta.get("n_discovery", "?"), meta.get("n_funnel", "?"),
                    meta.get("n_epitope", "?"), meta.get("n_survivors", "?")]]))

    # -- Grafik
    chart = outdir / "_chart.png"
    if _chart(peptides, chart):
        el.append(Spacer(1, 8))
        el.append(Image(str(chart), width=15 * cm, height=6.6 * cm))

    # -- En iyi adaylar
    el.append(Paragraph("2. En iyi aday peptitler (CDS kaynağı + lokus ile)", h2))
    rows = [["#", "Peptit", "Tip", "Skor", "CDS / kaynak", "Gen", "Lokus", "Allel"]]
    for i, p in enumerate(peptides[:15], 1):
        rows.append([i, p.seq, p.kind, f"{p.candidacy:.3f}",
                     str(p.parent)[:22],
                     str(p.metrics.get("gene") or "—"),
                     str(p.metrics.get("locus_tag") or "—"),
                     str(p.metrics.get("best_allele", "—"))])
    el.append(tbl(rows, widths=[0.7*cm, 3*cm, 1.3*cm, 1.2*cm, 3.6*cm, 1.7*cm, 1.9*cm, 2.2*cm]))

    # -- MHC anchor/cep motifi (yorum)
    mhci = [p for p in peptides if p.kind == "MHC-I" and p.metrics.get("anchor_residues")][:10]
    if mhci:
        el.append(Paragraph("2b. MHC yarığı anchor/cep motifi (yorum — sıralamayı etkilemez)", h2))
        el.append(Paragraph("Peptidin anchor kalıntıları (P2, C-terminal PΩ) ve o allelin "
                            "NetMHCpan taramasından AMPİRİK çıkarılan cep tercihi. Bağlanma uyumu "
                            "%rank'ta zaten puanlanır; bu tablo yalnız yorum içindir.", small))
        arows = [["Peptit", "Allel", "Anchorlar", "Allel cep tercihi", "Uyum"]]
        for p in mhci:
            anch = ";".join(f"{k}={v}" for k, v in (p.metrics.get("anchor_residues") or {}).items())
            motif = "  ".join(f"{k}:{','.join(v)}" for k, v in
                              (p.metrics.get("allele_anchor_motif") or {}).items()) or "—"
            arows.append([p.seq, str(p.metrics.get("best_allele", "—")), anch, motif,
                          str(p.metrics.get("anchor_match", "—"))])
        el.append(tbl(arows, widths=[2.6*cm, 2.8*cm, 2.6*cm, 5.2*cm, 1.2*cm]))

    # -- Kullanılan eşikler
    el.append(Paragraph("3. Kullanılan eşikler (tekrarlanabilirlik)", h2))
    trows = [["Adım", "Araç", "Parametre", "Değer", "Tip"]]
    for r in meta.get("thresholds", []):
        trows.append([r["step"], r["tool"], r["param"], f"{r['value']} {r['unit']}",
                      "sert" if r["hard_filter"] else "skor"])
    el.append(tbl(trows, widths=[2.6*cm, 3.4*cm, 3.2*cm, 3.4*cm, 1.6*cm]))

    # -- Yöntemler / araçlar (özet tablo)
    el.append(Paragraph("4. Kullanılan yöntemler ve araçlar", h2))
    refs = meta.get("citations") or citations.for_report()
    mrows = [["Adım", "Araç", "Ref"]] + [[r["step"], r["tool"], f"[{i}]"]
                                          for i, r in enumerate(refs, 1)]
    el.append(tbl(mrows, widths=[5*cm, 5*cm, 1.5*cm]))
    el.append(Spacer(1, 6))
    el.append(Paragraph("GPU gerektiren yapısal adımlar (AlphaFold peptit-MHC, moleküler dinamik) "
                        "bu koşuda ertelenmiştir (deferred).", small))

    # -- Popülasyon kapsamı (IEDB)
    popcov = meta.get("population_coverage") or {}
    if popcov:
        el.append(Paragraph("5. Popülasyon kapsamı (IEDB HLA frekansları)", h2))
        el.append(Paragraph("Aday epitop setinin, bir bireyin en az bir epitop-bağlayan "
                            "allele sahip olma olasılığı (%). Gerçek frekans verisi yalnız "
                            "insan HLA için mevcuttur.", small))
        areas = popcov.get("areas", [])
        for hname, he in popcov.get("hosts", {}).items():
            el.append(Spacer(1, 4))
            if he.get("note"):
                el.append(Paragraph(f"<b>{he.get('label', hname)}:</b> {he['note']}", small))
                continue
            prows = [["Sınıf"] + areas]
            for cls, lbl in (("mhc_i", "MHC-I"), ("mhc_ii", "MHC-II")):
                cov = he.get(cls)
                if not cov:
                    continue
                prows.append([lbl] + [f"{cov.get(a, {}).get('coverage', '—')}%"
                                      if a in cov else "—" for a in areas])
            if len(prows) > 1:
                el.append(Paragraph(f"<b>{he.get('label', hname)}</b>", body))
                el.append(tbl(prows))

    # -- Referanslar (tam atıflar)
    el.append(Paragraph("6. Referanslar", h2))
    for i, r in enumerate(refs, 1):
        doi = r["doi"]
        link = doi if doi.startswith("http") else f"https://doi.org/{doi}"
        el.append(Paragraph(f'[{i}] {r["citation"]} <font color="#1565c0">{link}</font>', small))
        el.append(Spacer(1, 2))

    doc.build(el)
    chart.unlink(missing_ok=True)
    return out
