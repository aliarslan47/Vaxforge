"""MANUSCRIPT(.md|_TR.md) -> .docx (tablo + figür gömülü).

Kullanım: python md_to_docx.py <girdi.md> <cikti.docx>
Biçim kuralı (kullanıcı): collapsible/Word-Heading stili YOK; bölüm etiketleri düz Normal
paragrafta bold + bir punto büyük (gövde 12pt, etiket 13pt, başlık 14pt). Markdown tabloları
gerçek docx tablosuna, ![](figures/..) resimleri gömülü görsele çevrilir. Em-dash zaten yok.
"""
import re
import sys
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

SRC = sys.argv[1] if len(sys.argv) > 1 else "MANUSCRIPT.md"
OUT = sys.argv[2] if len(sys.argv) > 2 else "MANUSCRIPT.docx"

doc = Document()
normal = doc.styles["Normal"]
normal.font.name = "Calibri"
normal.font.size = Pt(12)

_tok = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*)")
_img = re.compile(r"^!\[[^\]]*\]\(([^)]+)\)\s*$")

def clean(text):
    return text.replace("**", "").replace("*", "")

def add_runs(par, text, base_bold=False, base_italic=False, base_size=12):
    for part in _tok.split(text):
        if not part:
            continue
        bold, italic, inner = base_bold, base_italic, part
        if part.startswith("**") and part.endswith("**"):
            bold, inner = True, part[2:-2]
        elif part.startswith("*") and part.endswith("*"):
            italic, inner = True, part[1:-1]
        r = par.add_run(inner)
        r.bold = bold
        r.italic = italic
        r.font.size = Pt(base_size)

lines = open(SRC, encoding="utf-8").read().split("\n")
title_done = False
i = 0
n = len(lines)
while i < n:
    line = lines[i].rstrip()
    if line == "" or line == "---":
        i += 1
        continue

    # tablo bloğu (ardışık | ... | satırları)
    if line.startswith("|") and i + 1 < n and re.match(r"^\|[\s:|-]+\|", lines[i + 1]):
        header = [c.strip() for c in line.strip("|").split("|")]
        rows = []
        i += 2  # başlık + ayraç
        while i < n and lines[i].strip().startswith("|"):
            rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
            i += 1
        ncol = len(header)
        t = doc.add_table(rows=1, cols=ncol)
        t.style = "Light Grid Accent 1"
        for j, h in enumerate(header):
            cell = t.rows[0].cells[j]
            cell.paragraphs[0].text = ""
            add_runs(cell.paragraphs[0], clean(h), base_bold=True, base_size=11)
        for row in rows:
            cells = t.add_row().cells
            for j in range(ncol):
                txt = clean(row[j]) if j < len(row) else ""
                cells[j].paragraphs[0].text = ""
                bold = txt.startswith("Total") or txt.startswith("Toplam") or txt.startswith("VaxForge 14")
                add_runs(cells[j].paragraphs[0], txt, base_bold=bold, base_size=11)
        doc.add_paragraph()
        continue

    # gömülü figür
    m = _img.match(line)
    if m:
        path = m.group(1)
        try:
            doc.add_picture(path, width=Inches(6.2))
            doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        except Exception as e:
            p = doc.add_paragraph()
            add_runs(p, "[figür yüklenemedi: %s (%s)]" % (path, e))
        i += 1
        continue

    # tamamen **...** satır = başlık / bölüm etiketi
    if re.fullmatch(r"\*\*.+\*\*", line):
        inner = line[2:-2]
        p = doc.add_paragraph()
        if not title_done:
            add_runs(p, inner, base_bold=True, base_size=14)
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            title_done = True
        else:
            add_runs(p, inner, base_bold=True, base_size=13)
            p.paragraph_format.space_before = Pt(10)
        i += 1
        continue

    if line.startswith("- "):
        p = doc.add_paragraph(style="List Bullet")
        add_runs(p, line[2:])
        i += 1
        continue

    p = doc.add_paragraph()
    add_runs(p, line)
    i += 1

doc.save(OUT)
print("yazıldı:", OUT, "| paragraf:", len(doc.paragraphs), "| tablo:", len(doc.tables))
