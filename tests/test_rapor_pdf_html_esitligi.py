"""PDF ↔ HTML bölüm EŞİTLİĞİ (2026-11-23).

Bu repoda belgelenmiş dürüstlük kusurlarından biri "raporun SUSMASI"dır
(bkz. DURUM.md, on birinci varyant): bir bölüm HTML raporda vardı, PDF'te
hiç yoktu, dolayısıyla PDF'i okuyan kişi ölçümün yapıldığını bile görmüyordu.
Beşinci turda üç bölüm eklendi; 2026-11-23'te denetlendiğinde DÖRT bölüm hâlâ
yalnız HTML'deydi:

    glikozilasyon      → PDF'te hiç yok
    konservasyon       → PDF'te yalnız 'sınırlamalar' içinde tek satır
    moleküler mimikri  → PDF'te yalnız bir sayı satırı
    tanılama           → PDF'te hiç yok

Mimikri otoimmünite riskidir, konservasyon ise adaylık puanını besler; ikisinin
de tablosunun PDF'te bulunmaması bir ölçümün sessizce kaybolmasıdır.

Bu test o dördünü PDF'te ARAR. Amaç yalnız bugünü düzeltmek değil: bundan sonra
HTML'e bir bölüm eklenip PDF'e eklenmezse, asimetri burada yakalanır.
"""

import pytest

from vaxforge import report, report_pdf
from vaxforge.models import Peptide


def _pep(seq="TARFYRIFADIAQPP"):
    p = Peptide(seq=seq, parent="prot1", kind="MHC-II")
    p.passed, p.candidacy = True, 0.8
    p.notes = ["güçlü"]
    p.metrics["mhc_rank"] = 0.02
    p.metrics["self_mimicry"] = True
    return p


# Dört bölümün de VERİSİ dolu olduğu meta: bölüm ancak verisi varken basılır,
# bu yüzden testin "yok" demesi için verinin var olması şart.
_META = {
    "lang": "tr",
    "timestamp": "2026-11-23T00:00:00",
    "input": "test.gb",
    "input_path": "test.gb",
    "profile": "bacteria",
    "organism": None,
    "hosts": [],
    "n_input": 1, "n_discovery": 1, "n_funnel": 1, "n_epitope": 1, "n_survivors": 1,
    "conservation": {"computed": True, "n_strains": 5, "n_conserved": 3,
                     "min_percent": 80, "method": "GERÇEK (ortolog hizalama)"},
    "mimicry": {"available": True, "k": 9, "n_candidates": 25, "n_exact": 2,
                "n_near": 4, "hosts": ["human", "bovine"]},
    "glycosylation": {"available": True, "n_candidates": 25, "n_bcell": 11,
                      "n_bcell_masked": 3, "n_masked": 7},
    "diagnostics": [{"stage": "structure", "error": "HTTPError",
                     "message": "ESMAtlas 500", "where": "structure.py:88",
                     "expected": False}],
    "diagnostics_summary": {"unexpected": 1, "stages": ["structure"]},
}

# (bölüm anahtarı, PDF metninde aranacak parça, HTML üreten fonksiyon)
_BOLUMLER = [
    ("konservasyon", "Suşlar-arası konservasyon", report._conservation_section),
    ("mimikri", "Moleküler mimikri", report._mimicry_section),
    ("glikozilasyon", "Glikozilasyon", report._glyco_section),
    ("tanılama", "Tanılama", report._diagnostics_section),
]


@pytest.fixture(scope="module")
def pdf_metni(tmp_path_factory):
    pypdf = __import__("pypdf")
    out = report_pdf.build(tmp_path_factory.mktemp("rap"), [_pep()], dict(_META))
    return "\n".join((p.extract_text() or "") for p in pypdf.PdfReader(str(out)).pages)


@pytest.mark.parametrize("ad,baslik,html_fn", _BOLUMLER)
def test_html_bolumu_pdfte_de_var(ad, baslik, html_fn, pdf_metni):
    html = html_fn(dict(_META))
    assert html, f"{ad}: HTML bölümü boş döndü — test metası eksik"
    assert baslik in pdf_metni, (
        f"{ad} bölümü HTML raporda var, PDF'te YOK. Raporun susması da bir "
        f"dürüstlük kusurudur: iki çıktı biçimi eşit olmalı."
    )


def test_pdfte_mimikri_sayilari_tablo_olarak_var(pdf_metni):
    """Eskiden yalnız 'n adayda self-eşleşme' diye tek satır vardı; tam ve
    kısmi mimik SAYILARI ayrı ayrı görünmüyordu."""
    assert "Tam mimik" in pdf_metni
    assert "Kısmi mimik" in pdf_metni


def test_pdfte_konservasyon_olculmediginde_yuzde_yuz_denmez(tmp_path):
    """Konservasyon hesaplanmadıysa PDF 'Hesaplanmadı' yazar; %100 ima etmez."""
    pypdf = __import__("pypdf")
    meta = dict(_META)
    meta["conservation"] = {"computed": False}
    meta["conservation_auto_strains"] = {"note": "ortolog sayısı 2'nin altında"}
    out = report_pdf.build(tmp_path, [_pep()], meta)
    txt = "\n".join((p.extract_text() or "") for p in pypdf.PdfReader(str(out)).pages)
    assert "Hesaplanmadı" in txt
    assert "ortolog sayısı 2'nin altında" in txt


def test_pdfte_tanilama_beklenmeyen_hatayi_adiyla_yazar(pdf_metni):
    """'Araç kurulu değil' ile 'araç koşarken çöktü' ayrımı PDF'te de olmalı."""
    assert "HTTPError" in pdf_metni
    assert "structure" in pdf_metni


def test_konservasyon_olculemedi_sayisi_iki_raporda_da_gorunur(tmp_path):
    """csum['n_unmeasured'] run.json'da VARDI ama ne HTML ne PDF basıyordu:
    yüzdesi olmayan aday "korunmamış" değil, ÖLÇÜLEMEMİŞtir. 2026-11-22'deki
    gerçek koşuda 4 aday ölçülemedi ve bu hiçbir raporda görünmüyordu."""
    pypdf = __import__("pypdf")
    meta = dict(_META)
    meta["conservation"] = {"computed": True, "n_strains": 5, "n_conserved": 69,
                            "min_percent": 70.0, "n_unmeasured": 4,
                            "min_orthologs": 2, "n_identical_strains": 1,
                            "method": "GERÇEK (diamond ortolog)"}
    html = report._conservation_section(meta)
    assert "Ölçülemedi" in html and ">4<" in html
    out = report_pdf.build(tmp_path, [_pep()], meta)
    txt = "\n".join((p.extract_text() or "") for p in pypdf.PdfReader(str(out)).pages)
    assert "Ölçülemedi" in txt
    assert "Girdiyle aynı suş" in txt


def test_konservasyon_ozetinde_yontem_etiketi_var():
    """Rapor 'Yöntem: —' basıyordu: özet dict'te method anahtarı yoktu."""
    from vaxforge import conservation
    assert conservation.METHOD.startswith("GERÇEK")
