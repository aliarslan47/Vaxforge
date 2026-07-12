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

    # -- En iyi adaylar (hücreler Paragraph ile sarılır → taşma yok, uzun metin kaydırılır)
    el.append(Paragraph("2. En iyi aday peptitler (CDS kaynağı + lokus ile)", h2))
    cell = ParagraphStyle("cell", parent=body, fontSize=6.5, leading=8, wordWrap="CJK")

    def P(s):
        return Paragraph(str(s).replace("&", "&amp;").replace("<", "&lt;"), cell)

    rows = [["#", "Peptit", "Tip", "Skor", "CDS / kaynak", "Gen", "Lokus", "Allel"]]
    for i, p in enumerate(peptides[:15], 1):
        rows.append([str(i), P(p.seq), p.kind, f"{p.candidacy:.3f}",
                     P(str(p.parent)[:22]), P(p.metrics.get("gene") or "—"),
                     P(p.metrics.get("locus_tag") or "—"),
                     P(p.metrics.get("best_allele", "—"))])
    el.append(tbl(rows, widths=[0.7*cm, 3*cm, 1.2*cm, 1.1*cm, 3.2*cm, 1.6*cm, 2.6*cm, 3.0*cm]))

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

    # -- Aday epitoplar TİPE GÖRE sıralı tablolar (literatür stili: CTL/HTL/B)
    from . import evaluate
    tt = evaluate.type_tables(peptides, meta, top_n=15)
    total = {k: sum(1 for p in peptides if p.kind == k) for k in ("MHC-I", "MHC-II", "B")}
    el.append(Paragraph("2c. Aday epitoplar — tipe göre sıralı (CTL / HTL / B-hücre)", h2))
    el.append(Paragraph(
        "Epitoplar tipe göre ayrı tablolarda, her tip kendi içinde bağlanma gücüne göre sıralı "
        "(T-hücre: %rank artan; B-hücre: BepiPred azalan). ★ + yeşil satır = tüm zorunlu ölçütleri "
        "geçen final seçilen epitop (antijenik + güçlü bağlanma [+ HTL'de IFN-γ]). ✓=geçti, "
        "📖=IEDB literatürde, Kons.=suş verisi yok. Tam araç dökümü <b>candidates_full.xlsx</b>'te.", small))

    def _yn(b):
        return "✓" if b else "✗"

    def _n(v, nd=2):
        return "—" if v is None else f"{v:.{nd}f}"

    conf = {
        "MHC-I": ("🟦 CTL — MHC-I / CD8+ T-hücre",
                  ["★", "Epitop", "Kaynak", "%rank", "Allel", "Antij", "Al", "Tk", "İmmü", "İşle", "lit"],
                  [0.5*cm, 2.7*cm, 2.5*cm, 1.3*cm, 2.2*cm, 1.2*cm, 0.8*cm, 0.8*cm, 1.3*cm, 1.2*cm, 0.9*cm]),
        "MHC-II": ("🟨 HTL — MHC-II / CD4+ T-hücre",
                   ["★", "Epitop", "Kaynak", "%rank", "Allel", "Antij", "Al", "Tk", "IFN", "lit"],
                   [0.5*cm, 3.3*cm, 2.6*cm, 1.3*cm, 2.6*cm, 1.2*cm, 0.9*cm, 0.9*cm, 1.0*cm, 0.9*cm]),
        "B": ("🟩 B-hücre — antikor",
              ["★", "Epitop", "Kaynak", "Uz", "BepiPred", "Antij", "Al", "Tk", "lit"],
              [0.5*cm, 3.0*cm, 2.8*cm, 1.0*cm, 1.8*cm, 1.3*cm, 0.9*cm, 0.9*cm, 0.9*cm]),
    }
    for kind in ("MHC-I", "MHC-II", "B"):
        rows = tt.get(kind, [])
        if not rows:
            continue
        title, header, widths = conf[kind]
        el.append(Spacer(1, 6))
        el.append(Paragraph(f"{title} — {len(rows)}/{total[kind]} gösteriliyor", body))
        drows = [header]
        star_rows = []
        for ri, r in enumerate(rows, 1):
            star = "★" if r["star"] else ""
            lit = "📖" if r["iedb"] else "–"
            if r["star"]:
                star_rows.append(ri)
            if kind == "B":
                drows.append([star, r["epitope"], r["source"], str(r["length"]),
                              _n(r["bcell"]), _n(r["antigenicity"]),
                              _yn(r["allergen_ok"]), _yn(r["toxic_ok"]), lit])
            elif kind == "MHC-I":
                drows.append([star, r["epitope"], r["source"], _n(r["rank"]), r["allele"],
                              _n(r["antigenicity"]), _yn(r["allergen_ok"]), _yn(r["toxic_ok"]),
                              _n(r["immunogenicity"]), _n(r["processing"]), lit])
            else:
                drows.append([star, r["epitope"], r["source"], _n(r["rank"]), r["allele"],
                              _n(r["antigenicity"]), _yn(r["allergen_ok"]), _yn(r["toxic_ok"]),
                              _yn(r["ifn_ok"]), lit])
        t = tbl(drows, widths=widths)
        for ri in star_rows:
            t.setStyle(TableStyle([("BACKGROUND", (0, ri), (-1, ri), colors.HexColor("#e7f7ee"))]))
        el.append(t)

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

    # -- IEDB literatür/bilinen-epitop taraması + validasyon recall'ü
    im = meta.get("iedb_match")
    if im:
        el.append(Paragraph("5. IEDB literatür / bilinen-epitop taraması", h2))
        if not im.get("available"):
            el.append(Paragraph(im.get("note", "IEDB taraması yapılamadı."), small))
        else:
            el.append(Paragraph(
                f"Kaynak: {im.get('source','')}. Eşleşen aday: "
                f"<b>{im.get('n_matched',0)}/{im.get('n_candidates','?')}</b>. Bir adayın "
                "deneysel doğrulanmış bir epitopla örtüşmesi güçlü pozitif kontrol sinyalidir. "
                "Bu adım adaylık puanını değiştirmez (salt yorumlama).", small))
            mrows = [["#", "Peptit", "Tip", "Eşleşme", "IEDB epitobu", "Organizma", "PMID"]]
            for i, p in enumerate(peptides, 1):
                ie = p.metrics.get("iedb") or {}
                if ie.get("matched") is not True:
                    continue
                pmids = ie.get("pmids") or []
                ref = ", ".join(pmids[:3]) or (f"IEDB {ie.get('iedb_id')}" if ie.get("iedb_id") else "—")
                mrows.append([i, p.seq, p.kind, ie.get("match_type", ""),
                              ie.get("epitope_seq", ""),
                              (ie.get("organisms") or ["—"])[0], ref])
            if len(mrows) > 1:
                el.append(Spacer(1, 4))
                el.append(tbl(mrows, widths=[0.7*cm, 2.8*cm, 1.2*cm, 3.0*cm, 3.0*cm, 3.2*cm, 2.6*cm]))
            else:
                el.append(Paragraph("<i>Hiçbir aday bilinen IEDB epitobuyla eşleşmedi.</i>", small))
            bm = im.get("benchmark") or {}
            if bm and bm.get("recall") is not None:
                el.append(Spacer(1, 5))
                el.append(Paragraph("Validasyon — bilinen epitop recall'ü", body))
                el.append(Paragraph(
                    "Bu organizma için IEDB'deki deneysel doğrulanmış lineer epitoplar 'ground "
                    "truth' alınır; pipeline tahminleriyle örtüşme üzerinden ölçülür.", small))
                brows = [["Bilinen epitop", "Yakalanan", "Recall", "Tahmin", "Eşleşen aday", "Oran"],
                         [bm["n_known"], bm["n_known_hit"], f"{round(bm['recall']*100,1)}%",
                          bm["n_pred"], bm["n_pred_matched"],
                          f"{round((bm['precision_like'] or 0)*100,1)}%"]]
                el.append(tbl(brows))

    # -- Popülasyon kapsamı (IEDB)
    popcov = meta.get("population_coverage") or {}
    if popcov:
        el.append(Paragraph("6. Popülasyon kapsamı (IEDB HLA frekansları)", h2))
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
    el.append(Paragraph("7. Referanslar", h2))
    for i, r in enumerate(refs, 1):
        doi = r["doi"]
        link = doi if doi.startswith("http") else f"https://doi.org/{doi}"
        el.append(Paragraph(f'[{i}] {r["citation"]} <font color="#1565c0">{link}</font>', small))
        el.append(Spacer(1, 2))

    doc.build(el)
    chart.unlink(missing_ok=True)
    return out
