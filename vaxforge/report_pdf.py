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
from reportlab.platypus import (Image, PageBreak, Paragraph,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

from . import citations, correlate, mhc_motif
from .i18n import t

_FONT = "DejaVu"
_FONT_B = "DejaVu-Bold"

#: Sayfa kenar boşlukları (SimpleDocTemplate ile BİREBİR aynı olmak zorunda:
#: tablo genişlikleri bu çerçeveye göre ölçülür).
_MARGIN = 1.6 * cm
_FRAME_W = A4[0] - 2 * _MARGIN

#: Hücre iç boşluğu (pt). reportlab varsayılanı 6+6 pt = 0.42 cm idi; bu,
#: ölçülen metin genişliğine eklenince 'BoLA-DRB3_01201' gibi BÖLÜNMEMESİ
#: gereken simgeleri satır ortasından ikiye ayırıyordu (raporda allel adı
#: '...0120' + '1' olarak iki satıra düşüyordu).
_CELL_PAD = 2.2

#: Bir sütun tek başına çerçevenin en çok bu kadarını alabilir.
_MAX_COL_FRACTION = 0.45

#: Metin bu karakterlerden SONRA satır sonuna girebilir — yani bir simgenin
#: "bölünemez" en uzun parçası hesaplanırken bunlar ayırıcı sayılır.
_SOFT_BREAK = " /,;|()"


def _register_font() -> bool:
    base = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    try:
        pdfmetrics.registerFont(TTFont(_FONT, str(base / "DejaVuSans.ttf")))
        pdfmetrics.registerFont(TTFont(_FONT_B, str(base / "DejaVuSans-Bold.ttf")))
        return True
    except Exception:
        return False


def _chart(peptides, path: Path, lang: str = "tr") -> bool:
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
    ax.set_ylabel("Candidacy score" if lang=="en" else "Adaylık puanı")
    ax.set_title("Top candidate peptides (color = type)" if lang=="en" else "En iyi aday peptitler (renk = tip)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colmap.values()]
    ax.legend(handles, colmap.keys(), fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


# İzotip renkleri — frontend ISO_COLOR ile birebir aynı (görsel tutarlılık).
_ISO_COLOR = {
    "IgM": "#38bdf8", "IgG": "#34d399", "IgG1": "#34d399", "IgG2": "#a78bfa",
    "IgG2a": "#a78bfa", "IgG3": "#f59e0b", "IgA": "#f472b6", "IgE": "#fb7185",
    "IgY": "#34d399", "IgGa": "#34d399", "IgGb": "#a78bfa", "IgGT": "#f59e0b",
}


def _isotype_chart(curve: dict, path: Path, lang: str = "tr") -> bool:
    """Konağın simüle izotip-bazlı antikor yanıt eğrisini çizer (tepe=1 normalize).

    Veri run çıktısında zaten var (immunogenicity.per_host[i].curve); frontend
    kartıyla aynı grafiği PDF'e taşır.
    """
    days = (curve or {}).get("days") or []
    isotypes = (curve or {}).get("isotypes") or {}
    if len(days) < 2 or not isotypes:
        return False
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    for iso, ys in isotypes.items():
        if not ys:
            continue
        ax.plot(days[:len(ys)], ys, label=iso, color=_ISO_COLOR.get(iso, "#94a3b8"),
                linewidth=1.8, solid_joinstyle="round")
    ax.set_xlabel(t(lang, "pdf_immuno_day"), fontsize=8)
    ax.set_ylabel(t(lang, "pdf_immuno_relab"), fontsize=8)
    ax.set_title(t(lang, "pdf_immuno_curve"), fontsize=9)
    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc="upper right", ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def _pdf_ifn_cell(b, lang: str = "tr") -> str:
    """IFN-γ hücresi — ÜÇ DEĞERLİ. None = ÖLÇÜLMEDİ.

    '✗' basmak yasak: IFNepitope2 modeli yalnız insan/fare içindir, başka bir
    konakta metrik hiç hesaplanmaz ve '✗' okuyan kişi "IFN-γ indüklemiyor"
    sonucuna varırdı (ölçüm yokluğu olumsuz bulgu değildir).
    """
    if b is None:
        return "n/m" if lang == "en" else "ölç."
    return "✓" if b else "✗"


def _cell_text(c) -> str:
    """Hücrenin ÖLÇÜLECEK düz metni (Paragraph/Flowable ise metnini çıkarır)."""
    if isinstance(c, Paragraph):
        return getattr(c, "text", "") or ""
    return "" if c is None else str(c)


def _longest_token(s: str) -> str:
    """Satır sonuna bölünemeyen en uzun parça (bkz. _SOFT_BREAK)."""
    best, cur = "", ""
    for ch in s:
        cur += ch
        if ch in _SOFT_BREAK:
            if len(cur) > len(best):
                best = cur
            cur = ""
    return best if len(best) > len(cur) else cur


def _min_widths(data, header: bool, fsize: float, hsize: float) -> list[float]:
    """Sütun başına, hiçbir simgeyi ortadan BÖLMEYEN en küçük genişlik.

    Başlık satırı kalın ve bir punto büyük olduğu için ayrı ölçülür: raporda
    'iç içe girmiş' görünen başlıklar (ör. 'Bilinen epitop (benzersiz)' 3.6 cm,
    sütun 2.5 cm) tam olarak bu ölçümün hiç yapılmamasından doğuyordu.
    """
    nc = max(len(r) for r in data)
    mins = [0.0] * nc
    for ri, row in enumerate(data):
        bold = header and ri == 0
        fn, fs = (_FONT_B, hsize) if bold else (_FONT, fsize)
        for ci in range(nc):
            txt = _cell_text(row[ci] if ci < len(row) else "")
            tok = _longest_token(txt)
            try:
                w = pdfmetrics.stringWidth(tok, fn, fs)
            except Exception:
                w = fs * 0.6 * len(tok)
            mins[ci] = max(mins[ci], w + 2 * _CELL_PAD)
    return mins


def _fit_widths(data, widths, header: bool, fsize: float, hsize: float) -> list[float]:
    """Sütun genişliklerini çerçeveye ve içeriğe göre düzelt.

    Üç güvence (yeni raporların DAİMA düzgün çıkması için):
      1. Hiçbir sütun, bölünemez en uzun simgesinden dar kalmaz → başlık ya da
         hücre metni komşu sütunun üstüne TAŞMAZ.
      2. Toplam genişlik çerçeveyi (kenar boşlukları hariç sayfa) tam doldurur →
         ne sağa taşar ne de yarım sayfa boş kalır.
      3. Tek sütun çerçevenin %45'inden fazlasını yemez.
    """
    nc = max(len(r) for r in data)
    mins = _min_widths(data, header, fsize, hsize)
    cap = _FRAME_W * _MAX_COL_FRACTION
    w = [float(x) for x in (widths or [])][:nc]
    w += [_FRAME_W / nc] * (nc - len(w))
    w = [max(wi, min(mi, cap)) for wi, mi in zip(w, mins)]
    tot = sum(w) or 1.0
    return [wi * _FRAME_W / tot for wi in w]


def build(outdir: Path, peptides, meta: dict) -> Path:
    outdir = Path(outdir)
    lang = meta.get("lang", "tr")
    _register_font()
    out = outdir / "report.pdf"
    styles = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=styles["Normal"], fontName=_FONT, fontSize=9, leading=12)
    h1 = ParagraphStyle("h1", parent=styles["Heading1"], fontName=_FONT_B, fontSize=16, textColor=colors.HexColor("#0b6b4f"))
    h2 = ParagraphStyle("h2", parent=styles["Heading2"], fontName=_FONT_B, fontSize=11, spaceBefore=10)
    small = ParagraphStyle("small", parent=body, fontSize=7, textColor=colors.grey)
    mono = ParagraphStyle("mono", parent=body, fontName="Courier", fontSize=6, leading=7)

    doc = SimpleDocTemplate(str(out), pagesize=A4, topMargin=1.5 * cm,
                            bottomMargin=1.5 * cm,
                            leftMargin=_MARGIN, rightMargin=_MARGIN)
    el = []

    # Tablo hücre stilleri: BAŞLIK da dahil her hücre Paragraph'a sarılır.
    # Eskiden başlıklar çıplak str idi; reportlab çıplak str'i sarmaz, hücre
    # sınırından TAŞIRIR → komşu sütunun başlığıyla üst üste biniyordu.
    _TBL_FS, _TBL_HS = 7.0, 7.0
    # wordWrap="CJK" HER KARAKTERDE bölebilir: uzun peptit/DNA dizisi için
    # gerekli, ama boşluklu metinde sayıyı ortadan ikiye ayırıyordu (raporda
    # PMID '11163665' -> '1116366' + '5'). Bu yüzden iki stil: boşluk içeren
    # hücre KELİME sınırında, boşluksuz hücre her yerde bölünür.
    tcell = ParagraphStyle("tcell", parent=body, fontSize=_TBL_FS, leading=_TBL_FS + 1.6,
                           wordWrap="CJK")
    tcellw = ParagraphStyle("tcellw", parent=tcell, wordWrap=None)
    thead = ParagraphStyle("thead", parent=tcellw, fontName=_FONT_B,
                           fontSize=_TBL_HS, leading=_TBL_HS + 1.6,
                           textColor=colors.white)

    def _esc(s: str) -> str:
        return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))

    def tbl(data, widths=None, header=True):
        data = [list(r) for r in data]
        nc = max(len(r) for r in data)
        for r in data:
            r += [""] * (nc - len(r))
        widths = _fit_widths(data, widths, header, _TBL_FS, _TBL_HS)
        wrapped = []
        for ri, row in enumerate(data):
            cells = []
            for c in row:
                if hasattr(c, "wrapOn"):
                    cells.append(c)
                    continue
                txt = "" if c is None else str(c)
                st = (thead if (header and ri == 0)
                      else (tcellw if " " in txt.strip() else tcell))
                cells.append(Paragraph(_esc(txt), st))
            wrapped.append(cells)
        tb = Table(wrapped, colWidths=widths, repeatRows=1 if header else 0)
        style = [("FONTNAME", (0, 0), (-1, -1), _FONT), ("FONTSIZE", (0, 0), (-1, -1), _TBL_FS),
                 ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cccccc")),
                 ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                 ("LEFTPADDING", (0, 0), (-1, -1), _CELL_PAD),
                 ("RIGHTPADDING", (0, 0), (-1, -1), _CELL_PAD),
                 ("TOPPADDING", (0, 0), (-1, -1), 1.6),
                 ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
                 ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f6f6f6")])]
        if header:
            style += [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0b6b4f")),
                      ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                      ("FONTNAME", (0, 0), (-1, 0), _FONT_B)]
        tb.setStyle(TableStyle(style))
        return tb

    # Bölüm başlıkları: HER BAŞLIK KENDİ SAYFASINDAN BAŞLAR (kullanıcı kuralı).
    # Numara otomatik artar — eski elle yazılmış '2b/2c/2d' etiketleri bir bölüm
    # atlandığında yanlış/atlamalı numara veriyordu.
    _sec = {"n": 0}

    def sect(title: str, first: bool = False):
        _sec["n"] += 1
        if not first:
            el.append(PageBreak())
        el.append(Paragraph(f"{_sec['n']}. {title}", h2))

    def sub(title: str):
        """Alt başlık — kendi sayfasını açmaz, ama tablosundan ayrı düşmez."""
        el.append(Spacer(1, 6))
        el.append(Paragraph(f"<b>{title}</b>", body))

    # -- Başlık
    el.append(Paragraph(t(lang,"rep_title"), h1))
    hosts = ", ".join(h["label"] for h in meta.get("hosts", [])) or "—"
    el.append(Paragraph(f"{t(lang,'rep_generated')}: {meta.get('timestamp','')} &nbsp;·&nbsp; {t(lang,'rep_input')}: "
                        f"<b>{meta.get('input','')}</b> &nbsp;·&nbsp; {t(lang,'rep_profile')}: "
                        f"<b>{meta.get('profile','')}</b> &nbsp;·&nbsp; {t(lang,'rep_hosts')}: <b>{hosts}</b>", body))
    el.append(Spacer(1, 6))
    el.append(Paragraph(t(lang,"rep_disclaimer"), small))

    # -- İMMÜN KORELAT UYARILARI (ZORUNLU — özetten ÖNCE okunmalı)
    # Organizma, profilden değil girdiden çözülür: profil ('bacteria') organizmayı
    # söylemez, korelat ise organizmaya bağlıdır.
    _corr_org = (meta.get("organism")
                 or correlate.detect_from_input(meta.get("input_path"),
                                                extra=str(meta.get("input", ""))))
    _caveats = correlate.caveats(_corr_org)
    if _caveats:
        sect(t(lang, "corr_title"))
        el.append(Paragraph(t(lang, "corr_intro"), small))
        if any(c.severity == "critical" for c in _caveats):
            el.append(Spacer(1, 4))
            el.append(Paragraph(f"<b>{t(lang, 'corr_critical')}</b>", small))
        for c in _caveats:
            el.append(Spacer(1, 6))
            _sev = {"critical": "#b00020", "major": "#9a6700"}.get(c.severity, "#444")
            el.append(Paragraph(
                f'<font color="{_sev}"><b>[{c.severity.upper()}] {c.title(lang)}</b></font>',
                small))
            el.append(Paragraph(c.text(lang), small))
            if c.blocks_claim:
                el.append(Paragraph(
                    f"<i>{t(lang, 'corr_blocks')} \u201c{c.blocks_claim}\u201d</i>", small))
            for _cit in c.citations:
                el.append(Paragraph(f"&nbsp;&nbsp;• {_cit}", small))

    # -- Özet
    sect(t(lang, "rep_summary"))
    el.append(tbl([[t(lang,"rep_sum_input"), t(lang,"rep_sum_discovery"), t(lang,"rep_sum_funnel"), t(lang,"rep_sum_epitope"), t(lang,"rep_sum_survivors")],
                   [meta.get("n_input", "?"), meta.get("n_discovery", "?"), meta.get("n_funnel", "?"),
                    meta.get("n_epitope", "?"), meta.get("n_survivors", "?")]]))

    # -- Grafik
    chart = outdir / "_chart.png"
    _imm_charts: list[Path] = []   # 2e izotip-eğrisi geçici görselleri (sonda silinir)
    if _chart(peptides, chart, lang):
        el.append(Spacer(1, 8))
        el.append(Image(str(chart), width=15 * cm, height=6.6 * cm))

    # -- En iyi adaylar (hücreler Paragraph ile sarılır → taşma yok, uzun metin kaydırılır)
    sect(t(lang, "pdf_top"))

    def P(s):
        """Hücre metni. Sarma/kaçırma artık tbl() içinde, hücrenin içeriğine
        göre seçilen stille yapılıyor; burada yalnız metne çevirilir."""
        return "" if s is None else str(s)

    rows = [["#", t(lang,"col_epitope"), t(lang,"col_type"), "Score" if lang=="en" else "Skor", "CDS", t(lang,"col_source"), "Locus" if lang=="en" else "Lokus", t(lang,"col_allele")]]
    for i, p in enumerate(peptides[:15], 1):
        rows.append([str(i), P(p.seq), p.kind, f"{p.candidacy:.3f}",
                     P(str(p.parent)[:22]), P(p.metrics.get("gene") or "—"),
                     P(p.metrics.get("locus_tag") or "—"),
                     P(p.metrics.get("best_allele", "—"))])
    el.append(tbl(rows, widths=[0.7*cm, 3*cm, 1.2*cm, 1.1*cm, 3.2*cm, 1.6*cm, 2.6*cm, 3.0*cm]))

    # -- MHC anchor/cep motifi (yorum)
    # MHC-II de dâhil: anchor'lar her iki sınıfta da allele özgü ölçülüyor.
    _anchor_src = [p for p in peptides
                   if p.kind in ("MHC-I", "MHC-II") and p.metrics.get("anchor_residues")]
    mhci = _anchor_src[:10]
    if mhci:
        sect(t(lang, "pdf_anchor_title"))
        el.append(Paragraph(t(lang,"pdf_anchor_text"), small))
        def _anchor_full(m):
            # 'N/N' (tüm anchor'lar cebe uyuyor) = önemli → yıldız
            try:
                a, b = str(m).split("/")
                return int(a) == int(b) and int(b) > 0
            except Exception:
                return False

        arows = [["★", t(lang, "col_epitope"), t(lang, "col_allele"),
                  t(lang, "col_anchor_pos"), t(lang, "col_anchors"),
                  t(lang, "col_pocket"), t(lang, "col_match2")]]
        astar = []
        for p in mhci:
            anch = ";".join(f"{k}={v}" for k, v in (p.metrics.get("anchor_residues") or {}).items())
            motif = "  ".join(f"{k}:{','.join(v)}" for k, v in
                              (p.metrics.get("allele_anchor_motif") or {}).items()) or "—"
            am = str(p.metrics.get("anchor_match", "—"))
            full = _anchor_full(am)
            if full:
                astar.append(len(arows))
            arows.append(["★" if full else "", P(p.seq), P(str(p.metrics.get("best_allele", "—"))),
                          P(str(p.metrics.get("anchor_positions", "—"))),
                          P(anch), P(motif), am])
        atb = tbl(arows, widths=[0.5*cm, 2.3*cm, 2.3*cm, 1.6*cm, 2.1*cm, 4.3*cm, 1.0*cm])
        for ri in astar:
            atb.setStyle(TableStyle([("BACKGROUND", (0, ri), (-1, ri), colors.HexColor("#e7f7ee"))]))
        el.append(atb)
        # TABLONUN ALTI: anchor'lar allele özgü olduğundan, hangi TÜRÜN hangi
        # alleline göre bakıldığı yazılmak ZORUNDADIR; yoksa tablo yorumlanamaz.
        basis = mhc_motif.basis_lines(_anchor_src)
        if basis:
            el.append(Paragraph(f"<b>{t(lang, 'pdf_anchor_basis')}</b>", small))
            for b in basis:
                el.append(Paragraph(f"• {b}", small))
            el.append(Paragraph(t(lang, "pdf_anchor_basis_note"), small))

    # -- Aday epitoplar TİPE GÖRE sıralı tablolar (literatür stili: CTL/HTL/B)
    from . import evaluate
    tt = evaluate.type_tables(peptides, meta, top_n=15)
    total = {k: sum(1 for p in peptides if p.kind == k) for k in ("MHC-I", "MHC-II", "B")}
    sect(t(lang, "rep_bytype"))
    el.append(Paragraph(
        t(lang,"tt_intro_pdf"), small))

    def _yn(b):
        return _pdf_ifn_cell(b, lang)

    def _haz(state, value=None):
        """Tehlike sütunu: durumun kendisi (bkz. report.py hazard()).

        ✓/✗ kaldırıldı: ✓ "alerjen değil" demekti ama sütun başlığı tehlikenin
        adını taşıdığı için ters okunuyordu.
        """
        # PDF sütunu dar: 'alan dışı' uzun yazısı yerine '?' + skor. Anlamı
        # tablo altındaki dipnotta (tt_intro_pdf) açıklanır.
        if state == "ood":
            return "—" if value is None else f"? {value:.2f}"
        txt = t(lang, f"haz_{state}")
        return txt if value is None else f"{txt} {value:.2f}"

    def _n(v, nd=2):
        return "—" if v is None else f"{v:.{nd}f}"

    # Not: emoji (🟦/📖) DejaVu'da kutu çıkar → PDF'te metin başlık + ✓ kullanılır.
    # Hücreler Paragraph ile sarılır (taşma yok); kolonlar geniş (sayfa yeri var).
    _lit = "✓"
    conf = {
        "MHC-I": (t(lang, "tt_ctl"),
                  ["★", t(lang, "col_epitope"), t(lang, "col_source"), "%rank", t(lang, "col_allele"),
                   t(lang, "col_antig"), t(lang, "col_allergen"), t(lang, "col_toxic"),
                   t(lang, "col_immuno"), t(lang, "col_proc"), "Lit"],
                  [0.6*cm, 2.5*cm, 2.0*cm, 1.2*cm, 2.1*cm, 1.3*cm, 1.8*cm, 1.8*cm, 1.4*cm, 1.3*cm, 0.9*cm]),
        "MHC-II": (t(lang, "tt_htl"),
                   ["★", t(lang, "col_epitope"), t(lang, "col_source"), "%rank", t(lang, "col_allele"),
                    t(lang, "col_antig"), t(lang, "col_allergen"), t(lang, "col_toxic"), "IFN-γ", "Lit"],
                   [0.6*cm, 3.0*cm, 2.2*cm, 1.3*cm, 2.5*cm, 1.3*cm, 1.8*cm, 1.8*cm, 1.2*cm, 0.9*cm]),
        "B": (t(lang, "tt_bcell"),
              ["★", t(lang, "col_epitope"), t(lang, "col_source"), t(lang, "col_len"), "BepiPred",
               t(lang, "col_antig"), t(lang, "col_allergen"), t(lang, "col_toxic"), "Lit"],
              [0.6*cm, 3.2*cm, 2.4*cm, 1.3*cm, 1.8*cm, 1.4*cm, 1.8*cm, 1.8*cm, 0.9*cm]),
    }
    for kind in ("MHC-I", "MHC-II", "B"):
        rows = tt.get(kind, [])
        if not rows:
            continue
        title, header, widths = conf[kind]
        el.append(Spacer(1, 6))
        el.append(Paragraph(f"<b>{title}</b> — {len(rows)}/{total[kind]} {t(lang,'disc_short')}", body))
        drows = [header]
        star_rows = []
        for ri, r in enumerate(rows, 1):
            star = "★" if r["star"] else ""
            lit = _lit if r["iedb"] else "–"
            if r["star"]:
                star_rows.append(ri)
            epi, src, all = P(r["epitope"]), P(r["source"]), P(r.get("allele", ""))
            if kind == "B":
                drows.append([star, epi, src, str(r["length"]),
                              _n(r["bcell"]), _n(r["antigenicity"]),
                              _haz(r["allergen_state"]), _haz(r["toxic_state"], r.get("toxicity")), lit])
            elif kind == "MHC-I":
                drows.append([star, epi, src, _n(r["rank"]), all,
                              _n(r["antigenicity"]), _haz(r["allergen_state"]),
                              _haz(r["toxic_state"], r.get("toxicity")),
                              _n(r["immunogenicity"]), _n(r["processing"]), lit])
            else:
                drows.append([star, epi, src, _n(r["rank"]), all,
                              _n(r["antigenicity"]), _haz(r["allergen_state"]),
                              _haz(r["toxic_state"], r.get("toxicity")),
                              _yn(r["ifn_ok"]), lit])
        tb = tbl(drows, widths=widths)
        for ri in star_rows:
            tb.setStyle(TableStyle([("BACKGROUND", (0, ri), (-1, ri), colors.HexColor("#e7f7ee"))]))
        el.append(tb)
        if kind == "MHC-II" and any(r["ifn_ok"] is None for r in rows):
            el.append(Spacer(1, 3))
            el.append(Paragraph("⚠ " + t(lang, "ifn_note"), small))

    # refs bir kez hesaplanır — hem Yöntemler tablosu hem Referanslar kullanır.
    # SIRA (Ali, 2026-11-22): REFERANSLAR EN SONDA. Biyolojik kapsam/sınırlamalar
    # eşiklerden ÖNCE; IEDB taraması tip-tabloların HEMEN ARDINDAN.
    refs = meta.get("citations") or citations.for_report()

    # -- IEDB literatür/bilinen-epitop taraması + validasyon recall'ü
    im = meta.get("iedb_match")
    if im:
        sect(t(lang, "iedb_title"))
        if not im.get("available"):
            el.append(Paragraph(im.get("note", t(lang,"iedb_unavail")), small))
        else:
            el.append(Paragraph(
                t(lang,"iedb_intro").format(src=im.get('source',''), n_val=im.get('n_validated',0),
                                            n_weak=im.get('n_weak',0), tot=im.get('n_candidates','?')), small))
            el.append(Paragraph(t(lang,"iedb_caveat"), small))
            mrows = [["#", t(lang, "col_epitope"), t(lang, "col_type"), t(lang, "col_strength"),
                      t(lang, "col_iedb_epi"), t(lang, "col_organism"), t(lang, "col_assay"), "Ref"]]
            matched = [(i, p, p.metrics.get("iedb") or {})
                       for i, p in enumerate(peptides, 1)
                       if (p.metrics.get("iedb") or {}).get("matched") is True]
            matched.sort(key=lambda t3: (0 if t3[2].get("evidence") == "validated" else 1,
                                         -(t3[2].get("overlap") or 0)))
            for i, p, ie in matched:
                pmids = ie.get("pmids") or []
                ref = ", ".join(pmids[:2])
                if ie.get("iedb_id"):
                    ref = (ref + f", IEDB {ie.get('iedb_id')}") if ref else f"IEDB {ie.get('iedb_id')}"
                org = (ie.get("organisms") or ["—"])[0].split(" (")[0]   # uzun suş parantezini at
                assay = "; ".join((ie.get("assay_evidence") or [])[:1]) or "—"
                ev = ie.get("evidence")
                strength = (t(lang, "strength_validated") if ev == "validated"
                            else t(lang, "strength_weak")) + f" ({ie.get('overlap','')})"
                mrows.append([str(i), P(p.seq), p.kind, P(strength),
                              P(ie.get("epitope_seq", "")), P(org), P(assay), P(ref or "—")])
            if len(mrows) > 1:
                el.append(Spacer(1, 4))
                el.append(tbl(mrows, widths=[0.6*cm, 2.4*cm, 0.9*cm, 1.9*cm, 2.5*cm, 2.8*cm, 2.2*cm, 2.2*cm]))
            else:
                el.append(Paragraph(f"<i>{t(lang,'iedb_none')}</i>", small))
            bm = im.get("benchmark") or {}
            if bm and bm.get("recall") is not None:
                el.append(Spacer(1, 5))
                el.append(Paragraph(t(lang,"iedb_val_title"), body))
                el.append(Paragraph(
                    t(lang,"iedb_val_text").format(k=bm.get("k",8)), small))
                rex = bm.get("recall_exact")
                rex_s = f"{round(rex*100,1)}%" if rex is not None else "—"
                brows = [[t(lang,"col_known"), t(lang,"col_captured"), t(lang,"col_recall"), t(lang,"col_recall_exact"), t(lang,"col_pred"), t(lang,"col_pred_match"), t(lang,"col_matchrate")],
                         [bm["n_known"], bm["n_known_hit"], f"{round(bm['recall']*100,1)}%", rex_s,
                          bm["n_pred"], bm["n_pred_matched"],
                          f"{round((bm['precision_like'] or 0)*100,1)}%"]]
                el.append(tbl(brows))
                if bm.get("top_ref_share") and bm.get("top_ref_share") >= 0.5:
                    el.append(Paragraph(t(lang,"iedb_ref_dominance").format(
                        share=round(bm['top_ref_share']*100), pmid=bm.get('top_ref') or '—'), small))

    # -- KÜRATÖRLÜ DOĞRULAMA SETİ (organizmaya özgü; IEDB boş olduğunda)
    # IEDB bazı organizmalarda pratikte boştur (Leptospira: ~16 kayıt) →
    # IEDB tabanlı recall o organizmada ANLAMSIZDIR. Payda burada açıkça
    # yazılır ve KONAĞA GÖRE ayrılır.
    _cv = meta.get("curated_validation") or {}
    if _cv:
        sect(t(lang, "cv_title"))
        el.append(Paragraph(f"<b>{t(lang, 'cv_why')}</b> {_cv.get('why', '')}", small))
        el.append(Paragraph(f"<b>{t(lang, 'cv_scope')}</b> {_cv.get('scope_note', '')}",
                            small))
        el.append(Spacer(1, 4))
        el.append(tbl([[t(lang, "cv_n_total"), t(lang, "cv_n_pos"),
                        t(lang, "cv_n_neg"), t(lang, "cv_matched"),
                        t(lang, "cv_version")],
                       [_cv.get("n_total", 0), _cv.get("n_positive", 0),
                        _cv.get("n_negative", 0),
                        f"{_cv.get('n_matched', 0)}/{_cv.get('n_candidates', 0)}",
                        f"v{_cv.get('version')} ({_cv.get('compiled')})"]]))
        # Payda konağa göre
        _hc = _cv.get("host_coverage") or {}
        if _hc:
            el.append(Spacer(1, 4))
            el.append(Paragraph(f"<b>{t(lang, 'cv_by_host')}</b>", small))
            _hrows = [[t(lang, "cv_host"), t(lang, "cv_n_total"),
                       t(lang, "cv_n_pos"), t(lang, "cv_n_neg")]]
            for _h, _d in sorted(_hc.items(), key=lambda kv: -kv[1].get("n_total", 0)):
                _hrows.append([_h, _d.get("n_total", 0), _d.get("n_positive", 0),
                               _d.get("n_negative", 0)])
            el.append(tbl(_hrows, widths=[4.0*cm, 3.0*cm, 3.0*cm, 3.0*cm]))
        # Eşleşmeler
        _hits = _cv.get("hits") or []
        if _hits:
            el.append(Spacer(1, 4))
            _mrows = [[t(lang, "cv_candidate"), t(lang, "cv_epitope"),
                       "Protein", t(lang, "cv_match_type"), t(lang, "cv_label"),
                       t(lang, "cv_host"), "PMID"]]
            for _hit in _hits[:15]:
                _mrows.append([P(_hit["candidate"]), P(_hit["name"]),
                               P(_hit["protein"]), P(_hit["match_type"]),
                               P(_hit["label"]), P(_hit["host_serum"]),
                               P(_hit["pmid"])])
            el.append(tbl(_mrows, widths=[2.6*cm, 2.2*cm, 2.4*cm, 2.0*cm,
                                          1.8*cm, 2.6*cm, 1.4*cm]))
        else:
            el.append(Paragraph(t(lang, "cv_no_match"), small))
        _miss = _cv.get("missing_sequences") or {}
        if _miss:
            el.append(Spacer(1, 4))
            el.append(Paragraph(f"<b>{t(lang, 'cv_missing')}</b>", small))
            for _v in _miss.values():
                el.append(Paragraph(f"&nbsp;&nbsp;• {_v}", small))
        el.append(Paragraph(f"PMID: {', '.join(_cv.get('pmids', []))}", small))

    # -- PROTEİN düzeyi etiket seti (Wunder 2021 bağışıklık korelatları).
    # Yukarıdaki küratörlü set PEPTİT düzeyinde; hattın asıl çıktısı protein
    # sıralaması. Ölçüp basmamak bu repoda bilinen bir kusur (on ikinci
    # varyant), bu yüzden HTML ile AYNI turda PDF'e de eklendi.
    _ps = meta.get("protein_label_set") or {}
    if _ps:
        is_en = lang == "en"
        sect("Protein-level label set (correlates of immunity)" if is_en
             else "Protein düzeyi etiket seti (bağışıklık korelatları)")
        el.append(Paragraph(
            ("⚠ <b>Readout:</b> " if is_en else "⚠ <b>Etiketin okuması:</b> ")
            + P(_ps.get("readout", "")), small))
        el.append(Spacer(1, 3))
        el.append(Paragraph(P(_ps.get("citation", "")), small))
        el.append(Spacer(1, 4))
        el.append(Paragraph(
            (f"Set: <b>{_ps.get('n_features','—')}</b> array features / "
             f"<b>{_ps.get('n_proteins','—')}</b> unique proteins · locus tags "
             f"found in this run: {_ps.get('n_run_locus_tags', 0)}" if is_en else
             f"Set: <b>{_ps.get('n_features','—')}</b> array özelliği / "
             f"<b>{_ps.get('n_proteins','—')}</b> benzersiz protein · bu koşuda "
             f"bulunan locus_tag: {_ps.get('n_run_locus_tags', 0)}"), small))
        if _ps.get("note"):
            el.append(Spacer(1, 4))
            el.append(Paragraph("<i>" + P(_ps["note"]) + "</i>", small))
        else:
            el.append(Spacer(1, 4))
            el.append(Paragraph(
                ("The denominator is HOST-SPECIFIC: missing a hamster-significant "
                 "target in a mouse scan is not a failure. 'In input' separates "
                 "'was never in the input' from 'dropped by the filters'." if is_en else
                 "Payda KONAĞA GÖRE ayrıdır: hamsterda anlamlı bir hedefi fare "
                 "taramasında kaçırmak başarısızlık değildir. 'Girdide' sütunu, "
                 "'girdide hiç yoktu' durumunu 'elemede düştü' durumundan ayırır."),
                small))
            _with_input = any(d.get("in_input") is not None
                              for d in (_ps.get("hosts") or {}).values())
            _ph = [("Host" if is_en else "Konak"),
                   ("Denominator" if is_en else "Payda"),
                   ("Recovered" if is_en else "Yakalanan"),
                   ("Recall %" if is_en else "Recall %"),
                   ("Not meas." if is_en else "Ölçülmemiş")]
            if _with_input:
                _ph += [("In input" if is_en else "Girdide"),
                        ("Kept" if is_en else "Geçen"), "%"]
            prows3 = [_ph]
            for _h, _d in (_ps.get("hosts") or {}).items():
                _row = [_h, str(_d.get("denominator", "—")),
                        str(_d.get("recovered", "—")),
                        "—" if _d.get("recall_percent") is None
                        else str(_d["recall_percent"]),
                        str(_d.get("unmeasured_in_host", "—"))]
                if _with_input:
                    _row += [str(_d.get("in_input", "—")),
                             str(_d.get("kept_of_in_input", "—")),
                             "—" if _d.get("kept_percent_of_in_input") is None
                             else str(_d["kept_percent_of_in_input"])]
                prows3.append(_row)
            el.append(Spacer(1, 4))
            el.append(tbl(prows3))
            if _ps.get("specificity_note"):
                el.append(Spacer(1, 3))
                el.append(Paragraph("<i>" + P(_ps["specificity_note"]) + "</i>", small))

    # -- Aşı Konstruktu (MEV)
    mv = meta.get("mev")
    if mv:
        p = mv["properties"]
        is_en = lang == "en"
        ag, al, tx = p["antigenicity"], p["allergen"], p["toxicity"]
        hs, ss, ds, so = (p["human_similarity"], p.get("secondary_structure", {}),
                          p.get("disorder", {}), p.get("solubility", {}))

        def _yn3(v, yes, no):
            return "—" if v is None else (yes if v else no)

        nb = mv["n_by_kind"]
        sect("Vaccine Construct (MEV)" if is_en else "Aşı Konstruktu (MEV)")
        el.append(Paragraph(
            (f"β-defensin adjuvant + EAAAK/AAY/GPGPG/KK linkers · "
             f"CTL={nb['MHC-I']}, HTL={nb['MHC-II']}, LBL={nb['B']} · "
             f"length: <b>{len(mv['sequence'])} aa</b>" if is_en
             else f"β-defensin adjuvan + EAAAK/AAY/GPGPG/KK linker · "
             f"CTL={nb['MHC-I']}, HTL={nb['MHC-II']}, LBL={nb['B']} · "
             f"uzunluk: <b>{len(mv['sequence'])} aa</b>"), small))
        el.append(Paragraph(mv["sequence"], mono))
        el.append(Spacer(1, 3))
        mrows = [
            ["Property" if is_en else "Özellik", "Value" if is_en else "Değer"],
            ["MW / pI", f"{p['mw_kda']} kDa / {p['pI']}"],
            [("Charged −/+ (D+E / R+K)" if is_en else "Yüklü −/+ (D+E / R+K)"),
             f"{p['neg_residues']} / {p['pos_residues']}"],
            ["Instability", f"{p['instability']} ({'Stable' if p['stable'] else 'Unstable'})"
             if is_en else f"{p['instability']} ({'Kararlı' if p['stable'] else 'Kararsız'})"],
            ["Aliphatic / GRAVY", f"{p['aliphatic_index']} / {p['gravy']}"],
            ["Antigenicity", f"{ag['score']} · {ag['method']}"],
            ["Allergenicity", _yn3(al["allergenic"], "Allergen", "Non-allergen") + f" · {al['method']}"],
            ["Toxicity", _yn3(tx["toxic"], "Toxic", "Non-toxic") + f" · {tx['method']}"],
            [("Human similarity" if is_en else "İnsan benzerliği"),
             ("No similarity" if is_en else "Benzerlik yok") if hs.get("available") and not hs.get("similar")
             else (f"{hs.get('best_pident')}%" if hs.get("available") else "—")],
            ["Solubility" if is_en else "Çözünürlük",
             f"{so.get('scaled_solubility')} ({_yn3(so.get('soluble'),'Soluble' if is_en else 'Çözünür','Low')})"
             if so.get("scaled_solubility") is not None else "—"],
            [("Disordered / α-Helix / β-Strand" if is_en else "Düzensiz / α-Heliks / β-Şerit"),
             f"{ds.get('percent_disordered')}% / {ss.get('helix_pct')}% / {ss.get('strand_pct')}%"
             if ds.get("percent_disordered") is not None else "—"],
        ]
        el.append(tbl(mrows, widths=[5.5*cm, 8.7*cm]))

        # Alerjenite ATFI — "Allergen" satırı tek başına YANILTICI: eşleşmenin
        # konstruktun NERESİNDEN doğduğunu söylemeden okuyan kişi epitopların
        # alerjen olduğu sonucuna varır. Oysa bu koşuda üç eşleşmenin üçü de
        # linker/epitop BİRLEŞME YERİNDE doğdu (epitopların hepsi temiz) ve
        # linker değiştirilerek giderilebilir. Bkz. mev.allergen_attribution.
        hits = al.get("matches") or []
        if hits:
            sub("Allergenicity attribution (FAO/WHO 6-mer)" if is_en
                else "Alerjenite atfı (FAO/WHO 6-mer)")
            reg_tr = {"junction": "birleşme yeri", "linker": "linker",
                      "adjuvant": "adjuvan", "epitope": "epitop", "other": "diğer"}
            hrows = [[("6-mer" if is_en else "6-mer"), ("Pos." if is_en else "Konum"),
                      ("Region" if is_en else "Bölge"), ("Parts" if is_en else "Parçalar"),
                      ("Matched allergen" if is_en else "Benzediği alerjen")]]
            for h in hits:
                reg = h.get("region", "other")
                hrows.append([P(h.get("kmer", "")), str(h.get("start", "")),
                              reg if is_en else reg_tr.get(reg, reg),
                              P(h.get("detail", "—") or "—"), P(h.get("source", "—") or "—")])
            el.append(tbl(hrows, widths=[2.0*cm, 1.4*cm, 2.4*cm, 4.2*cm, 4.2*cm]))
            by = al.get("by_region") or {}
            if al.get("from_epitope"):
                verdict = ("At least one match lies INSIDE an epitope: the candidate filter "
                           "missed it — drop that epitope from the construct." if is_en
                           else "En az bir eşleşme bir EPİTOBUN İÇİNDE: aday filtresi kaçırmış — "
                           "o epitop konstruktan çıkarılmalı.")
            elif al.get("fixable_by_linker"):
                verdict = ("All matches arise at linker/epitope junctions — no epitope is itself "
                           "allergenic; changing the linker (or the order of epitopes) removes them."
                           if is_en else
                           "Eşleşmelerin tümü linker/epitop birleşme yerinde doğuyor — hiçbir "
                           "epitop kendiliğinden alerjen DEĞİL; linker (ya da epitop sırası) "
                           "değiştirilerek giderilebilir.")
            else:
                verdict = ("Matches lie in the host-derived adjuvant; this is expected and is not "
                           "a property of the selected epitopes." if is_en else
                           "Eşleşmeler konak kaynaklı adjuvanda; bu BEKLENEN bir durumdur ve "
                           "seçilen epitopların bir özelliği değildir.")
            el.append(Spacer(1, 2))
            el.append(Paragraph(
                ("Regions: " if is_en else "Bölgeler: ")
                + ", ".join(f"{(k if is_en else reg_tr.get(k, k))}={v}" for k, v in by.items())
                + " · " + verdict, small))

        # Düşük karmaşıklık nedeniyle KONSTRUKTA KONMAYAN epitoplar — sessizce
        # düşürülmezler, gerekçesiyle yazılırlar (bkz. complexity.py).
        exc = mv.get("excluded_low_complexity") or []
        if exc:
            sub("⚠ " + t(lang, "lc_title"))
            el.append(Paragraph(t(lang, "lc_intro"), small))
            erows = [[t(lang, "col_epitope"), t(lang, "col_type"), t(lang, "col_source"),
                      "%rank", ("Score" if is_en else "Skor"),
                      ("Reason" if is_en else "Gerekçe")]]
            for e in exc[:20]:
                erows.append([P(e.get("seq", "")), e.get("kind", ""), P(e.get("source", "—")),
                              "—" if e.get("rank") is None else f"{e['rank']:.3f}",
                              "—" if e.get("candidacy") is None else f"{e['candidacy']:.3f}",
                              P(e.get("reason", ""))])
            el.append(tbl(erows, widths=[3.2*cm, 1.2*cm, 2.4*cm, 1.3*cm, 1.3*cm, 5.0*cm]))

    # -- İmmünojenisite (ImmForge — in silico eleme; DÜRÜST etiketli)
    imm = meta.get("immunogenicity")
    if imm and imm.get("per_host"):
        is_en = lang == "en"
        sect("Immunogenicity — in silico screen (ImmForge)" if is_en
             else "İmmünojenisite — in silico eleme (ImmForge)")
        el.append(Paragraph(
            ("Mechanistic engine PREDICTION — not in vivo efficacy, NOT validated. Only the kinetic "
             "shape was validated on peptide-vaccine data; the verdict/isotype/profile are operational "
             "(threshold WORKS≥60). Isotype identities come from literature steering (primary IgM→IgG "
             "subclass mixture); relative magnitudes are operational (depend on the Th1/Th2 balance, "
             "not validated per construct)." if is_en
             else "Mekanistik motor ÖNGÖRÜSÜDÜR — in vivo etkinlik DEĞİL, valide EDİLMEMİŞTİR. Yalnızca "
             "kinetik şekil, peptit-aşısı verisiyle doğrulanmıştır; karar/profil operasyoneldir "
             "(eşik: OLUR≥60). İzotip kimlikleri literatür yönlendirmesinden gelir (birincil IgM→IgG "
             "alt-sınıf karışımı); göreli büyüklükler operasyoneldir (Th1/Th2 dengesine bağlıdır, "
             "konstrukt başına valide edilmemiştir)."), small))
        el.append(Spacer(1, 3))
        vmap = {"OLUR": ("WORKS" if is_en else "OLUR"), "ZAYIF": ("WEAK" if is_en else "ZAYIF"),
                "OLMAZ": ("FAILS" if is_en else "OLMAZ")}
        rows = [[("Host" if is_en else "Konak"), ("Verdict" if is_en else "Karar"),
                 ("Score" if is_en else "Skor"), "MHC-I", "MHC-II",
                 ("Clearance" if is_en else "Temizlik"), ("Profile" if is_en else "Profil"), "Ig"]]
        for ph in imm["per_host"]:
            pr = ph.get("profile") or {}
            rows.append([
                P(ph.get("host_label", ph.get("host", ""))),
                vmap.get(ph.get("verdict", ""), ph.get("verdict", "")),
                f"{round(ph.get('score', 0))}/100", str(ph.get("n_mhci", "")), str(ph.get("n_mhcii", "")),
                f"{round(ph.get('cleared_fraction', 0) * 100)}%",
                P(f"{pr.get('polarization', '—')} ({round(pr.get('th1', 0) * 100)}/{round(pr.get('th2', 0) * 100)})"),
                pr.get("dominant_isotype") or "—",
            ])
        el.append(tbl(rows, widths=[3.0*cm, 1.8*cm, 1.6*cm, 1.3*cm, 1.3*cm, 1.9*cm, 2.5*cm, 1.4*cm]))
        r0 = (imm["per_host"][0].get("reasons") or [])
        if r0:
            el.append(Spacer(1, 2))
            el.append(Paragraph("· " + " · ".join(r0), small))
        # Simüle antikor yanıt eğrisi — her konak için (frontend kartıyla aynı grafik).
        for hi, ph in enumerate(imm["per_host"]):
            icurve = outdir / f"_immchart_{hi}.png"
            if _isotype_chart(ph.get("curve") or {}, icurve, lang):
                el.append(Spacer(1, 4))
                el.append(Paragraph(
                    f"<b>{ph.get('host_label', ph.get('host', ''))}</b>", small))
                el.append(Image(str(icurve), width=13 * cm, height=4.7 * cm))
                _imm_charts.append(icurve)

    # -- Yapısal değerlendirme (ESMFold/pLDDT) — PDF'te HİÇ YOKTU.
    # HTML raporda vardı, PDF'te yoktu: pLDDT 47.7 ("çok düşük — düzensiz
    # model") olan bir koşuda PDF, yapıdan türeyen hiçbir sonucun güvenini
    # göstermiyordu. Rapor eksikliği de bir dürüstlük kusurudur.
    st = meta.get("structure")
    if st:
        is_en = lang == "en"
        sect("Structural validation — ESMFold" if is_en
             else "Yapısal doğrulama — ESMFold")
        el.append(Paragraph(
            ("ESMFold PREDICTION via ESMAtlas API — not an experimental structure. "
             "pLDDT>70 confident, <50 disordered. Conformational-B is DiscoTope-LIKE "
             "(RSA+Parker), not the official DiscoTope." if is_en else
             "ESMAtlas API üzerinden ESMFold ÖNGÖRÜSÜ — deneysel yapı değil. pLDDT>70 "
             "güvenilir, <50 düzensiz. Konformasyonel-B, DiscoTope-BENZERİ (RSA+Parker); "
             "resmi DiscoTope değil."), small))
        if not st.get("available"):
            el.append(Paragraph(st.get("note") or ("Structure not predicted." if is_en
                                                   else "Yapı öngörülmedi."), small))
        else:
            surf = st.get("surface") or {}
            cb = st.get("conformational_b") or {}
            cl = st.get("plddt_confidence") or {}
            thr = st.get("plddt_threshold", 70.0)
            lbl = cl.get("label_en" if is_en else "label_tr", "")
            if cl and not cl.get("reliable", True):
                el.append(Spacer(1, 3))
                el.append(Paragraph(
                    f"<b>⚠ pLDDT {st.get('plddt_mean','—')} &lt; {thr} — {lbl}.</b> "
                    + ("Everything read OFF this model is low-confidence too: the surface % "
                       "(RSA) and the conformational B-cell patches below. This is EXPECTED "
                       "for a MEV — the construct is a linker-joined chimera, not a natural "
                       "folding unit — so a low pLDDT is not by itself a defect, but the "
                       "surface/conformational results cannot carry the weight of a "
                       "measurement." if is_en else
                       "Bu modelden OKUNAN her şey de düşük güvenlidir: aşağıdaki yüzey % "
                       "(RSA) ve konformasyonel B-hücre yamaları. Bu, MEV için BEKLENEN "
                       "durumdur — konstrükt linkerlarla birleştirilmiş yapay bir kimera, "
                       "doğal bir katlanma birimi değil; yani düşük pLDDT tek başına kusur "
                       "değil, ama yüzey/konformasyonel sonuçlar ölçüm ağırlığı taşımaz."),
                    small))
            tag = "" if cl.get("reliable", True) else (" (low confidence)" if is_en
                                                      else " (düşük güven)")
            srows = [[("Property" if is_en else "Özellik"), ("Value" if is_en else "Değer")],
                     [("Source" if is_en else "Kaynak"), P(st.get("source", "—"))],
                     [("Tool" if is_en else "Araç"), P(st.get("tool", "—"))],
                     [("Residues" if is_en else "Kalıntı"), str(st.get("n_res", "—"))],
                     ["pLDDT (mean/min)" if is_en else "pLDDT (ort/min)",
                      P(f"{st.get('plddt_mean','—')} / {st.get('plddt_min','—')}"
                        + (f" — {lbl}" if lbl else ""))],
                     [("Surface %" if is_en else "Yüzey %"),
                      f"{surf.get('percent_surface','—')}%{tag}"],
                     [("Conformational B patches" if is_en else "Konformasyonel B yaması"),
                      f"{cb.get('n_patches','—')}{tag}"]]
            el.append(tbl(srows, widths=[5.5*cm, 8.7*cm]))
            patches = cb.get("patches") or []
            if patches:
                sub("Surface B-cell patches" if is_en else "Yüzey B-hücre yamaları")
                prows = [[("Range" if is_en else "Aralık"), ("Len" if is_en else "Uzunluk"),
                          ("Peptide" if is_en else "Peptit"), ("Score" if is_en else "Skor"),
                          "pLDDT", ("Confidence" if is_en else "Güven")]]
                lvl_tr = {"very_high": "çok yüksek", "confident": "güvenilir",
                          "low": "düşük", "disordered": "düzensiz",
                          "unknown": "bilinmiyor"}
                for q in patches[:12]:
                    _pl = q.get("plddt_mean")
                    lv = str(q.get("plddt_confidence", "—"))
                    prows.append([f"{q['start']}–{q['end']}", str(q["length"]),
                                  P(q["peptide"]), str(q["mean_score"]),
                                  "—" if _pl is None else str(_pl),
                                  ("⚠ " if q.get("low_confidence") else "")
                                  + (lv if is_en else lvl_tr.get(lv, lv))])
                el.append(tbl(prows, widths=[1.8*cm, 1.4*cm, 4.4*cm, 1.6*cm, 1.6*cm, 3.0*cm]))

    # -- Suşlar-arası konservasyon (HTML'de TAM bölüm, PDF'te yalnız tek satır
    # 'sınırlamalar' içinde duruyordu). Konservasyon adaylık puanını besliyor ve
    # 2026-11'de iki kusuru düzeltildi (tek-suş %100 yanılgısı + profil kayması);
    # yeni dürüst değerler ("ölçülemedi", ortolog sayısı) PDF'e geçmiyordu.
    cons_ = meta.get("conservation") or {}
    if cons_:
        is_en = lang == "en"
        sect("Cross-strain conservation" if is_en else "Suşlar-arası konservasyon")
        if not cons_.get("computed"):
            auto = meta.get("conservation_auto_strains") or {}
            note = auto.get("note") or ("no strain data" if is_en else "suş verisi yok")
            el.append(Paragraph(
                ("<b>Not computed</b> — " if is_en else "<b>Hesaplanmadı</b> — ")
                + P(note), small))
            el.append(Spacer(1, 3))
            el.append(Paragraph(
                ("A conservation value is reported only when at least two orthologs "
                 "were actually compared. Where it could not be measured the field is "
                 "empty, never 100%." if is_en else
                 "Konservasyon değeri ancak EN AZ İKİ ortolog gerçekten "
                 "karşılaştırıldığında raporlanır. Ölçülemediği yerde alan boştur, "
                 "asla %100 değildir."), small))
        else:
            crows2 = [[("Property" if is_en else "Özellik"), ("Value" if is_en else "Değer")],
                      [("Strains" if is_en else "Suş sayısı"), str(cons_.get("n_strains", "—"))],
                      [("Conserved candidates" if is_en else "Korunmuş aday"),
                       str(cons_.get("n_conserved", "—"))],
                      [("Threshold" if is_en else "Eşik"),
                       f"≥%{cons_.get('min_percent', '—')}"],
                      # Ölçüm YAPILAMAYAN aday sayısı: yüzdesi olmayan aday
                      # "korunmamış" değildir, ÖLÇÜLEMEMİŞTİR.
                      [("Not measurable (fewer than min. orthologs)" if is_en
                        else "Ölçülemedi (en az ortolog sayısının altında)"),
                       str(cons_.get("n_unmeasured", "—"))],
                      [("Min. orthologs required" if is_en else "Gereken en az ortolog"),
                       str(cons_.get("min_orthologs", "—"))],
                      [("Strains identical to input (excluded)" if is_en
                        else "Girdiyle aynı suş (sayılmadı)"),
                       str(cons_.get("n_identical_strains", "—"))],
                      [("Method" if is_en else "Yöntem"),
                       P(cons_.get("method") or "—")]]
            if cons_.get("auto_fetched"):
                crows2.append([("Auto-fetched" if is_en else "Otomatik çekilen"),
                               P(str(cons_.get("auto_note", "—")))])
            el.append(tbl(crows2, widths=[5.5*cm, 8.7*cm]))

    # -- Moleküler mimikri (HTML'de TAM bölüm, PDF'te yalnız bir sayı satırı).
    # Bu bölüm OTOİMMÜNİTE riskidir: PDF'i okuyan kişi kaç adayın konak
    # proteiniyle tam k-mer eşleştiğini tablo olarak göremiyordu.
    mim_ = meta.get("mimicry") or {}
    if mim_.get("available"):
        is_en = lang == "en"
        sect("Molecular mimicry — host self-match" if is_en
             else "Moleküler mimikri — konak self-eşleşme")
        el.append(Paragraph(
            ("⚠ Exact k-mer match = strong autoimmunity risk; BLOSUM62 Hamming-1 "
             "(single conservative substitution) = partial cross-reactivity risk. "
             "SOFT — this flags candidates, it does not filter them." if is_en else
             "⚠ Tam k-mer eşleşme = güçlü otoimmünite riski; BLOSUM62 Hamming-1 "
             "(tek konservatif ikame) = kısmi çapraz-reaktivite riski. YUMUŞAK — "
             "bu ölçüt adayı işaretler, ELEMEZ."), small))
        el.append(Spacer(1, 4))
        mrows2 = [[("Property" if is_en else "Özellik"), ("Value" if is_en else "Değer")],
                  [("k-mer length" if is_en else "k-mer uzunluğu"), str(mim_.get("k", "—"))],
                  [("Candidates" if is_en else "Aday"), str(mim_.get("n_candidates", "—"))],
                  [("Exact mimics (hard)" if is_en else "Tam mimik (sert)"),
                   str(mim_.get("n_exact", "—"))],
                  [("Partial mimics (BLOSUM Hamming-1)" if is_en
                    else "Kısmi mimik (BLOSUM Hamming-1)"), str(mim_.get("n_near", "—"))],
                  [("Hosts" if is_en else "Konaklar"),
                   P(", ".join(mim_.get("hosts", [])) or "—")]]
        el.append(tbl(mrows2, widths=[5.5*cm, 8.7*cm]))

    # -- Glikozilasyon (HTML'de vardı, PDF'te HİÇ YOKTU)
    gly_ = meta.get("glycosylation") or {}
    if gly_.get("available"):
        is_en = lang == "en"
        sect("Glycosylation — glycan-masked epitopes" if is_en
             else "Glikozilasyon — glikan-maskeli epitoplar")
        el.append(Paragraph(
            ("⚠ N-X-S/T is an exact sequon motif, but actual glycan occupancy needs "
             "experimental proof; the O-glyco figure here is a Ser/Thr-rich HEURISTIC. "
             "PREDICTION — masked B-cell epitopes may be antibody-inaccessible in vivo "
             "(soft flag, nothing filtered)." if is_en else
             "⚠ N-X-S/T kesin bir sequon motifidir, ama gerçek glikan işgali deneysel "
             "olarak doğrulanmalıdır; buradaki O-gliko sayısı Ser/Thr-yoğun bir "
             "HEURİSTİKTİR. ÖNGÖRÜ — maskeli B-epitopları in vivo antikora "
             "erişilemeyebilir (yumuşak uyarı; aday elenmedi)."), small))
        el.append(Spacer(1, 4))
        grows = [[("Property" if is_en else "Özellik"), ("Value" if is_en else "Değer")],
                 [("Candidates scanned" if is_en else "Taranan aday"),
                  str(gly_.get("n_candidates", "—"))],
                 [("B-cell epitopes" if is_en else "B-hücre epitobu"),
                  str(gly_.get("n_bcell", "—"))],
                 [("B-cell glycan-masked" if is_en else "B-hücre glikan-maskeli"),
                  str(gly_.get("n_bcell_masked", "—"))],
                 [("All masked (N-sequon)" if is_en else "Tüm maskeli (N-sequon)"),
                  str(gly_.get("n_masked", "—"))]]
        el.append(tbl(grows, widths=[5.5*cm, 8.7*cm]))

    # -- In-siliko klonlama (HTML'de vardı, PDF'te yoktu)
    cl_ = meta.get("cloning")
    if cl_ and cl_.get("available"):
        is_en = lang == "en"
        sect("In-silico cloning" if is_en else "In-siliko klonlama")
        crows = [[("Property" if is_en else "Özellik"), ("Value" if is_en else "Değer")],
                 [("Host" if is_en else "Konak"), P(cl_.get("host", "—"))],
                 [("Method" if is_en else "Yöntem"), P(cl_.get("method", "—"))],
                 [("Insert" if is_en else "İnsert"),
                  f"{cl_.get('dna_len_bp','—')} bp / {cl_.get('protein_len','—')} aa"],
                 ["GC %", str(cl_.get("gc_percent", "—"))],
                 ["CAI", f"{cl_.get('cai','—')} (baseline {cl_.get('cai_baseline','—')})"],
                 [("Vector / tag" if is_en else "Vektör / etiket"),
                  P(f"{cl_.get('suggested_vector','—')} · {cl_.get('tag','—')}")],
                 [("Internal restriction sites" if is_en else "İç kesim bölgeleri"),
                  P("; ".join(f"{k}:{v}" for k, v in
                              (cl_.get("restriction_sites") or {}).items()) or "—")]]
        el.append(tbl(crows, widths=[5.5*cm, 8.7*cm]))
        if cl_.get("note"):
            el.append(Spacer(1, 3))
            el.append(Paragraph(cl_["note"], small))

    # -- Araç durumu (21 araç) — PDF'te HİÇ YOKTU. Bir aracın kurulu olmadığı
    # ya da doğrulama koşusunun düştüğü, PDF'i okuyan kişiye görünmüyordu.
    ts = meta.get("tool_status") or {}
    tsum = meta.get("tool_status_summary") or {}
    if ts:
        is_en = lang == "en"
        sect("Tool status" if is_en else "Araç durumu")
        el.append(Paragraph(
            (f"{tsum.get('available','?')}/{tsum.get('total','?')} tools present "
             f"({tsum.get('verified','?')} verified by a probe run, "
             f"{tsum.get('shallow','?')} file/import check only). A missing tool is "
             "separated into 'not installed' and 'present but its verification run "
             "failed' — the second can be a transient load error, and calling it "
             "'absent' would be just as misleading as claiming it ran." if is_en else
             f"{tsum.get('available','?')}/{tsum.get('total','?')} araç yerinde "
             f"({tsum.get('verified','?')}'ü sonda koşusuyla doğrulandı, "
             f"{tsum.get('shallow','?')} tanesi yalnız dosya/import denetimi). Eksiklik iki "
             "türe ayrılır: 'kurulu değil' ve 'dosya var ama doğrulama koşusu düştü' — "
             "ikincisi geçici yük hatası olabilir; ona 'araç yok' demek, çalıştığını "
             "iddia etmek kadar yanıltıcıdır."), small))
        el.append(Spacer(1, 3))
        tr_chk = {"verified": "sonda koşusu", "shallow": "dosya/import"}
        trows2 = [[("Tool" if is_en else "Araç"), ("State" if is_en else "Durum"),
                   ("Check" if is_en else "Denetim"), ("Note" if is_en else "Not")]]
        for key, v in ts.items():
            okk = bool(v.get("available"))
            fk = v.get("fail_kind")
            state = (("present" if is_en else "var") if okk else
                     {"missing": ("not installed" if is_en else "kurulu değil"),
                      "verify_failed": ("verify failed" if is_en else "doğrulama düştü"),
                      "probe_error": ("probe error" if is_en else "sonda hatası")}
                     .get(fk, ("missing" if is_en else "yok")))
            chk = v.get("check", "—")
            trows2.append([P(v.get("label_en" if is_en else "label", key)), state,
                           chk if is_en else tr_chk.get(chk, chk),
                           P(v.get("note") or "—")])
        el.append(tbl(trows2, widths=[4.6*cm, 2.4*cm, 2.2*cm, 5.6*cm]))
        # Kurulu AMA bu koşuda uygulanamayan araç: IFN-γ (konak modeli yok).
        if any(r["ifn_ok"] is None for r in (tt.get("MHC-II") or [])):
            el.append(Spacer(1, 3))
            el.append(Paragraph("⚠ " + t(lang, "ifn_note"), small))




    # -- Tanılama: yakalanan istisnalar (HTML'de vardı, PDF'te HİÇ YOKTU).
    # Hat aşamaları try/except ile sarılı; bir araç çökünce koşu DEVAM eder.
    # Bu bölüm olmadan "araç kurulu değil" ile "araç koşarken çöktü" PDF'te
    # ayırt edilemiyordu — araç durumu bölümünün tamamlayıcısı.
    _diag_items = meta.get("diagnostics") or []
    _diag_summ = meta.get("diagnostics_summary") or {}
    if _diag_items or _diag_summ:
        is_en = lang == "en"
        sect("Diagnostics — caught exceptions" if is_en else "Tanılama — yakalanan istisnalar")
        _n_unexp = int(_diag_summ.get("unexpected") or 0)
        if not _n_unexp:
            el.append(Paragraph("✓ " + t(lang, "diag_clean"), small))
        else:
            _stages = ", ".join(_diag_summ.get("stages") or [])
            el.append(Paragraph(f"⚠ <b>{t(lang, 'diag_unexpected')}</b>: "
                                f"{_n_unexp} ({_stages})", small))
            el.append(Spacer(1, 4))
            _dh = (("Stage", "Error", "Message", "Where") if is_en
                   else ("Aşama", "Hata", "Mesaj", "Yer"))
            drows = [[_dh[0], _dh[1], _dh[2], _dh[3]]]
            for d in _diag_items:
                if d.get("expected"):
                    continue
                drows.append([P(d.get("stage", "—")), P(d.get("error", "—")),
                              P((d.get("message") or "—")[:160]),
                              P(d.get("where") or "—")])
            if len(drows) > 1:
                el.append(tbl(drows, widths=[3.0*cm, 3.0*cm, 6.2*cm, 2.0*cm]))

    # -- Popülasyon kapsamı (IEDB)
    popcov = meta.get("population_coverage") or {}
    if popcov:
        sect(t(lang, "pop_title"))
        el.append(Paragraph(t(lang,"pop_text_pdf"), small))
        areas = popcov.get("areas", [])
        for hname, he in popcov.get("hosts", {}).items():
            el.append(Spacer(1, 4))
            if he.get("note"):
                el.append(Paragraph(f"<b>{he.get('label', hname)}:</b> {he['note']}", small))
                continue
            prows = [[t(lang,"col_type")] + areas]
            for cls, lbl in (("mhc_i", "MHC-I"), ("mhc_ii", "MHC-II")):
                cov = he.get(cls)
                if not cov:
                    continue
                prows.append([lbl] + [f"{cov.get(a, {}).get('coverage', '—')}%"
                                      if a in cov else "—" for a in areas])
            if len(prows) > 1:
                el.append(Paragraph(f"<b>{he.get('label', hname)}</b>", body))
                el.append(tbl(prows))
    # -- Biyolojik kapsam ve dürüst sınırlamalar (konformasyonel B / konservasyon / mimikri)
    sect(t(lang, "pdf_bio_limits_title"))
    el.append(Paragraph("• "+t(lang,"pdf_bcell_caveat"), small))
    cons = meta.get("conservation") or {}
    if cons.get("computed"):
        el.append(Paragraph("• "+t(lang,"pdf_cons_computed").format(
            n=cons.get("n_strains", 0), c=cons.get("n_conserved", 0),
            m=cons.get("min_percent", 80)), small))
    else:
        el.append(Paragraph("• "+t(lang,"pdf_cons_skipped"), small))
    n_mim = sum(1 for p in peptides if p.metrics.get("self_mimicry"))
    el.append(Paragraph("• "+t(lang,"pdf_mimicry_note").format(n=n_mim), small))

    # -- Kullanılan eşikler (tekrarlanabilirlik) — biyolojik sınırlamalardan SONRA
    sect(t(lang, "rep_thresholds"))
    trows = [[t(lang,"col_step"), t(lang,"col_tool"), t(lang,"col_param"), t(lang,"col_value"), t(lang,"col_type")]]
    for r in meta.get("thresholds", []):
        trows.append([P(r["step"]), P(r["tool"]), P(r["param"]), P(f"{r['value']} {r['unit']}"),
                      t(lang,"type_hard") if r["hard_filter"] else t(lang,"type_score")])
    el.append(tbl(trows, widths=[2.8*cm, 3.6*cm, 3.4*cm, 3.2*cm, 1.8*cm]))

    # Eşik tablosu CONFIG NİYETİNİ gösterir; profile göre bir adım hiç
    # uygulanmamış olabilir (ör. virüste lokalizasyon elemesi yapılmaz, ama
    # tabloda 'sert' yazar). Bu koşuda FİİLEN kullanılan yöntemler:
    fm = meta.get("funnel_methods") or {}
    if fm:
        is_en = lang == "en"
        el.append(Spacer(1, 4))
        el.append(Paragraph(
            ("The table above is the CONFIGURED intent. Depending on the pathogen "
             "profile a step may not have been applied at all in this run (for a "
             "virus, for instance, there is no localization filter even though the "
             "row says 'hard'). What was actually used in this run:" if is_en else
             "Yukarıdaki tablo YAPILANDIRILMIŞ niyeti gösterir. Patojen profiline "
             "göre bir adım bu koşuda hiç uygulanmamış olabilir (ör. virüste "
             "lokalizasyon elemesi yapılmaz, satırda 'sert' yazsa bile). Bu koşuda "
             "FİİLEN kullanılan yöntemler:"), small))
        frows = [[("Step" if is_en else "Adım"), ("Method actually used" if is_en
                                                  else "Bu koşuda kullanılan yöntem")]]
        for k, v in fm.items():
            frows.append([P(k), P(v)])
        el.append(tbl(frows, widths=[4.4*cm, 9.8*cm]))

    # -- Kullanılan yöntemler ve araçlar (özet tablo) — eşiklerden sonra, referanslardan ÖNCE
    sect(t(lang, "rep_methods"))
    mrows = [[t(lang,"col_step"), t(lang,"col_tool"), "Ref"]] + [[P(r["step"]), P(r["tool"]), f"[{i}]"]
                                          for i, r in enumerate(refs, 1)]
    el.append(tbl(mrows, widths=[5.5*cm, 8.1*cm, 1.2*cm]))
    el.append(Spacer(1, 6))
    el.append(Paragraph(t(lang,"pdf_methods_note"), small))
    # -- Referanslar (tam atıflar) — RAPORUN EN SONU (Ali kuralı 2026-11-22).
    # Yöntemler tablosundaki [n] numaraları buraya ileri atıf yapar; makale
    # düzeninde olağan olan budur.
    sect(t(lang, "rep_references"))
    for i, r in enumerate(refs, 1):
        doi = r["doi"]
        link = doi if doi.startswith("http") else f"https://doi.org/{doi}"
        # Araç bu makinede YOK ise atıf satırı bunu söyler (HTML raporla aynı
        # kural): yoksa kaynakça, koşmamış bir aracı koşmuş gibi gösterir.
        absent = ("" if r.get("used") is not False
                  else f' <font color="#a8601a"><b>[{t(lang,"tool_status_absent")}]</b></font>')
        el.append(Paragraph(f'[{i}] {r["citation"]}{absent} '
                            f'<font color="#1565c0">{link}</font>', small))
        el.append(Spacer(1, 2))
    doc.build(el)
    chart.unlink(missing_ok=True)
    for ic in _imm_charts:
        ic.unlink(missing_ok=True)
    return out
