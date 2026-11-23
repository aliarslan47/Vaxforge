"""Wunder 2021 protein düzeyi etiket seti — yükleme ve dürüstlük kuralları."""

import pytest

from vaxforge import lepto_targets as lt


def test_etiket_seti_yuklenir_ve_sayilar_makaleyle_tutar():
    assert lt.available()
    # Makale "154 protein hedefi" diyor; tabloda 154 SATIR var ama 152
    # benzersiz ORF: LIC10464 (ligB) dizide üç ayrı parça olarak basılmış.
    assert lt.n_features() == 154
    assert lt.n_proteins() == 152


def test_paydasi_konaga_gore_ayridir():
    """Tek payda kullanmak yasak: anlamlılık konağa göre değişiyor.

    Payda BENZERSİZ PROTEİN sayar, satır değil: hamsterda 133 satır anlamlı
    ama 131 protein (ligB'nin üç parçası da anlamlı çıkıyor). Satır saymak
    paydayı şişirir ve recall'ı yapay olarak düşürür."""
    d = {h: len(lt.denominator(h)) for h in lt.HOSTS}
    assert d["hamster"] == 131   # 133 satır
    assert d["mouse"] == 54      # 56 satır
    assert d["mouse_all"] == 13
    # Üçü de farklı — yani "Leptospira'da recall" diye tek sayı yoktur.
    assert len(set(d.values())) == 3


def test_konak_zorunlu_ve_bilinmeyen_konak_reddedilir():
    with pytest.raises(ValueError):
        lt.denominator("human")
    with pytest.raises(ValueError):
        lt.unmeasured("bovine")


def test_olculmemis_hucre_negatife_cevrilmez():
    """Fare sütunlarında 18 satırda ölçüm YOK. Bunlar 'anlamsız' sayılırsa
    payda şişer ve recall yapay olarak düşer."""
    assert len(lt.unmeasured("mouse")) == 18
    assert len(lt.unmeasured("mouse_all")) == 18
    assert len(lt.unmeasured("hamster")) == 0
    # None ile False ayrı: hamsterda 21 satır ÖLÇÜLDÜ ve anlamsız çıktı.
    sig = [t.significant_in("hamster") for t in lt.load()]
    assert sig.count(False) == 21
    assert sig.count(None) == 0


def test_locus_tag_yalniz_acik_kalipla_eslesir():
    assert lt.locus_tags("LIC10464 Ig-like repeat") == {"LIC10464"}
    assert lt.locus_tags("lic_10215 something") == {"LIC10215"}
    # Ürün adı benzerliği eşleşme SAYILMAZ: 'Hypothetical protein' 85 satırda
    # geçiyor, ada göre eşleştirme uydurma olurdu.
    assert lt.locus_tags("Hypothetical protein") == set()
    assert lt.locus_tags(None) == set()


class _Pr:
    def __init__(self, pid, desc="", locus=None):
        self.id, self.description = pid, desc
        self.annotations = {"locus_tag": locus} if locus else {}


def test_eslesme_yoksa_recall_sifir_denmez():
    """Sıfır eşleşme çoğu zaman 'kimlik biçimi tutmadı' demektir; bunu
    'recall %0' diye raporlamak adayı haksız yere kötü gösterir."""
    r = lt.match([_Pr("sp|Q9XYZ1|TEST", "Hypothetical protein")])
    assert "note" in r
    assert "recall %0 DEĞİL" in r["note"]
    assert r["hosts"] == {}


def test_eslesme_varsa_her_konak_icin_ayri_recall():
    hits = [_Pr(o) for o in lt.denominator("mouse_all")[:5]]
    r = lt.match(hits + [_Pr("LIC99999")])
    assert r["n_run_locus_tags"] == 6
    assert r["hosts"]["mouse_all"]["recovered"] == 5
    assert r["hosts"]["mouse_all"]["denominator"] == 13
    assert r["hosts"]["mouse_all"]["recall_percent"] == pytest.approx(38.5, abs=0.1)
    # Aynı eşleşme kümesi, farklı konak → farklı payda ve farklı recall.
    assert r["hosts"]["hamster"]["denominator"] == 131
    assert r["hosts"]["hamster"]["recall_percent"] != \
        r["hosts"]["mouse_all"]["recall_percent"]
    assert "özgüllük ÖLÇÜLEMEZ" in r["specificity_note"]


def test_etiketin_koruma_olmadigi_metinde_yazili():
    """Bu organizmada reaktivite ile koruma kopuk (Barbosa 2024: MCP güçlü IgG
    uyardı, korumadı, hastalığı HIZLANDIRDI). Etiket adı bunu saklamamalı."""
    assert "korelat" in lt.READOUT
    assert "koruma ölçümü değil" in lt.READOUT
    assert "eLife" in lt.CITATION and "64166" in lt.CITATION
    assert lt.METHOD.startswith("GERÇEK")


# ── Rapora GİRİYOR mu? (bu repoda "ölçüp basmama" bilinen bir kusur) ────────
def _meta_with_label_set():
    from vaxforge import lepto_targets as _lt
    hits = [_Pr(o) for o in _lt.denominator("mouse_all")[:4]]
    inp = set(_lt.denominator("mouse_all")[:8])
    return {
        "lang": "tr", "timestamp": "2026-11-23T00:00:00",
        "input": "lepto.fasta", "input_path": "lepto.fasta", "profile": "bacteria",
        "organism": "leptospira", "hosts": [],
        "n_input": 1, "n_discovery": 1, "n_funnel": 1, "n_epitope": 1,
        "n_survivors": 1,
        "protein_label_set": _lt.match(hits, input_tags=inp),
    }


def test_protein_etiket_seti_html_raporda_var():
    from vaxforge import report
    html = report._protein_label_set_section(_meta_with_label_set())
    assert "Protein düzeyi etiket seti" in html
    assert "korelat" in html            # koruma OLMADIĞI yazılı
    assert "64166" in html              # atıf var
    assert "hamster" in html and "mouse_all" in html


def test_protein_etiket_seti_pdfte_de_var(tmp_path):
    pypdf = __import__("pypdf")
    from vaxforge import report_pdf
    from vaxforge.models import Peptide
    p = Peptide(seq="TARFYRIFADIAQPP", parent="prot1", kind="MHC-II")
    p.passed, p.candidacy, p.notes = True, 0.8, ["güçlü"]
    out = report_pdf.build(tmp_path, [p], _meta_with_label_set())
    txt = "\n".join((pg.extract_text() or "") for pg in pypdf.PdfReader(str(out)).pages)
    assert "Protein düzeyi etiket seti" in txt
    assert "Girdide" in txt             # 'girdide yoktu' ≠ 'elemede düştü'
    assert "özgüllük" in txt            # özgüllük ölçülemez uyarısı
