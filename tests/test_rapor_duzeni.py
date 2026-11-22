"""2026-11-22 beşinci tur — PDF rapor düzeni ve ölçüm-yokluğu kusurları.

Ali'nin raporda (run_2026-11-22T20-34-48) bulduğu/istediği dört şey:
  1. Her bölüm başlığı KENDİ SAYFASINDAN başlar.
  2. Tablo başlıkları iç içe girmez (ölçüldü: eski PDF'te 'Kaynak antijen'
     başlığı komşu sütunun başlığının 13.7 pt üzerine biniyordu).
  3. Hiçbir simge satır ortasından bölünmez (eski PDF: allel
     'BoLA-DRB3_0120' + '1', PMID '1116366' + '5').
  4. IFN-γ ölçülmediğinde '✗' yazılmaz ve ★ ölçütü engellenmez.
Ayrıca: rapor koşu TEKRARLANMADAN yeniden üretilebilir (peptides.json).
"""

import json

from reportlab.lib.units import cm

from vaxforge import evaluate, peptide_io, report_pdf
from vaxforge.models import Peptide


# ── 1. Sütun genişliği / taşma ─────────────────────────────────────────────
def test_longest_token_yumusak_ayiricilari_sayar():
    # Ayırıcı, kendinden ÖNCEKİ parçaya ait sayılır (satır sonunda o karakter
    # de yer kaplar) — ölçüm bu yüzden kasıtlı olarak muhafazakâr.
    assert report_pdf._longest_token("10067674, 11163665") == "10067674,"
    assert report_pdf._longest_token("ab 1116366512") == "1116366512"
    assert report_pdf._longest_token("BoLA-DRB3_01201") == "BoLA-DRB3_01201"
    assert report_pdf._longest_token("Bilinen epitop (benzersiz)") == "benzersiz)"


def test_fit_widths_cerceveyi_tam_doldurur():
    data = [["a", "b", "c"], ["1", "2", "3"]]
    w = report_pdf._fit_widths(data, [1 * cm, 1 * cm, 1 * cm], True, 7, 7)
    assert abs(sum(w) - report_pdf._FRAME_W) < 0.01


def test_fit_widths_bolunemez_simgeyi_sigdirir():
    """ESKİ HATA: sütun 2.5 cm, 'BoLA-DRB3_01201' + hücre boşluğu 2.57 cm →
    allel adı ortadan ikiye bölünüyordu."""
    data = [["Allel"], ["BoLA-DRB3_01201"]]
    w = report_pdf._fit_widths(data, [2.5 * cm], True, 7, 7)
    need = report_pdf._min_widths(data, True, 7, 7)[0]
    assert w[0] >= need


def test_fit_widths_tek_sutun_sayfayi_yemez():
    uzun = "X" * 400
    data = [["a", "b"], [uzun, "k"]]
    w = report_pdf._fit_widths(data, None, True, 7, 7)
    # Üst sınır BÜYÜMEYİ sınırlar: 400 karakterlik tek simge sütunu şişirip
    # tabloyu sayfadan taşıramaz; toplam yine çerçeveye eşittir.
    assert w[0] <= report_pdf._FRAME_W * 0.5 + 0.01
    assert abs(sum(w) - report_pdf._FRAME_W) < 0.01


# ── 2. IFN-γ: ölçüm yokluğu olumsuz bulgu değil ────────────────────────────
def _mhcii(seq="TARFYRIFADIAQPP", ifn=None, rank=0.02):
    p = Peptide(seq=seq, parent="prot1", kind="MHC-II")
    p.passed, p.candidacy = True, 0.8
    p.notes = ["güçlü"]
    p.metrics = {"parent_antigenicity": 0.55, "pseudo_rank": rank,
                 "best_allele": "BoLA-DRB3_01201", "allergen": False,
                 "toxicity": 0.5}
    if ifn is not None:
        p.metrics["ifn_gamma_inducer"] = ifn
    return p


_META = {"thresholds": [
    {"step": "antigen_funnel", "tool": "antigenicity_vaxijen",
     "param": "threshold", "value": 0.4, "unit": "", "hard_filter": False},
    {"step": "survival", "tool": "toxicity", "param": "max_toxicity_score",
     "value": 0.6, "unit": "", "hard_filter": True}]}


def test_ifn_olculmediyse_ucuncu_durum():
    r = evaluate.type_tables([_mhcii(ifn=None)], _META)["MHC-II"][0]
    assert r["ifn_ok"] is None, "ölçülmemiş değer False olarak sunulamaz"
    assert report_pdf._pdf_ifn_cell(r["ifn_ok"], "tr") != "✗"


def test_ifn_olculmediyse_yildiz_engellenmez():
    """ESKİ HATA: ★ ölçütü `ifn_gamma_inducer is True` istiyordu; sığır
    konağında metrik hiç ölçülmediği için 209 HTL epitobunun hiçbiri ★
    alamadı (raporda HTL tablosunda tek yıldız yoktu)."""
    assert evaluate.type_tables([_mhcii(ifn=None)], _META)["MHC-II"][0]["star"] is True


def test_ifn_olculup_negatifse_yildiz_engellenir():
    assert evaluate.type_tables([_mhcii(ifn=False)], _META)["MHC-II"][0]["star"] is False
    assert evaluate.type_tables([_mhcii(ifn=True)], _META)["MHC-II"][0]["star"] is True


# ── 3. Her başlık kendi sayfasında (gerçek PDF üretilip okunur) ────────────
def _build_pdf(tmp_path):
    peps = [_mhcii(seq="TARFYRIFADIAQPP"), _mhcii(seq="PRPYVLVAGHEAAPF")]
    meta = dict(_META)
    meta.update({"lang": "tr", "timestamp": "2026-11-22T20:34:48",
                 "input": "test.gb", "profile": "virus", "hosts": [],
                 "n_input": 1, "n_discovery": 1, "n_funnel": 1,
                 "n_epitope": 2, "n_survivors": 2,
                 "tool_status": {"a": {"label": "Araç A", "label_en": "Tool A",
                                       "available": True, "check": "shallow",
                                       "note": ""}},
                 "tool_status_summary": {"total": 1, "available": 1,
                                         "verified": 0, "shallow": 1}})
    return report_pdf.build(tmp_path, peps, meta)


def test_her_baslik_kendi_sayfasinda(tmp_path):
    """Ali'nin kuralı: 'raporda her başlık kendi sayfasından başlayacak'."""
    pypdf = __import__("pypdf")
    out = _build_pdf(tmp_path)
    pages = pypdf.PdfReader(str(out)).pages
    import re
    seen = 0
    for pg in pages:
        lines = [ln.strip() for ln in (pg.extract_text() or "").splitlines()
                 if ln.strip()]
        for ln in lines[1:]:
            assert not re.match(r"^\d+\.\s+\S", ln), \
                f"başlık sayfanın ORTASINDA: {ln!r}"
        if lines and re.match(r"^\d+\.\s+\S", lines[0]):
            seen += 1
    assert seen >= 4, f"numaralı bölüm başlığı sayfa başında bulunamadı ({seen})"


def test_numaralandirma_atlamaz(tmp_path):
    """Eski kod '2b/2c/2d' etiketlerini ELLE yazıyordu: bir bölüm atlanınca
    numara atlıyor/yanlış çıkıyordu. Artık otomatik ve kesintisiz."""
    pypdf = __import__("pypdf")
    import re
    out = _build_pdf(tmp_path)
    txt = "\n".join((p.extract_text() or "") for p in pypdf.PdfReader(str(out)).pages)
    nums = [int(m) for m in re.findall(r"^(\d+)\.\s+\S", txt, re.M)]
    assert nums == list(range(1, len(nums) + 1)), nums


# ── 4. Peptit dökümü: rapor koşu tekrarlanmadan üretilebilir ───────────────
def test_peptides_json_gidis_donus(tmp_path):
    p = _mhcii()
    p.metrics["iedb"] = {"matched": True, "overlap": 9}
    p.methods["mhc_score"] = "GERÇEK (NetMHCIIpan)"
    out = peptide_io.dump([p], tmp_path / "peptides.json")
    assert json.loads(out.read_text())["version"] == peptide_io.VERSION
    back = peptide_io.load(out)[0]
    assert (back.seq, back.kind, back.parent) == (p.seq, p.kind, p.parent)
    assert back.metrics == p.metrics and back.methods == p.methods
    assert back.notes == p.notes and back.passed is True
