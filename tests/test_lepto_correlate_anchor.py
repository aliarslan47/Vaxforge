"""Leptospira doğrulama seti, korelat uyarıları ve allele özgü anchor testleri.

Hiçbiri ağ ya da NetMHCpan GEREKTİRMEZ: araç yoksa dürüst düşüş yolları test
edilir (anchor pozisyonu UYDURULMADIĞI doğrulanır).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from vaxforge import correlate, mhc_motif
from vaxforge import lepto_validation as lv

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


# --------------------------------------------------------------------------- #
# Doğrulama seti
# --------------------------------------------------------------------------- #
def test_koordinatlar_uniprot_dizisiyle_birebir():
    """coords alanları UniProt O34094'e karşı doğrulanmalı (dizi uydurma kontrolü)."""
    problems = lv.verify_coords()
    assert problems == [], "koordinat uyuşmazlığı: " + "; ".join(problems)


def test_referans_dizi_kaynaktan_okunuyor():
    seq = lv.lipl32_seq()
    assert len(seq) == 272, "UniProt O34094 272 aa olmalı"
    assert set(seq) <= AA20


def test_her_peptit_kanit_tasiyor():
    """Atıfsız/kaynaksız peptit sete girmemeli."""
    for p in lv.all_peptides():
        assert p.seq and set(p.seq.upper()) <= AA20, f"{p.name}: geçersiz dizi"
        assert p.pmid and p.pmid.isdigit(), f"{p.name}: PMID yok"
        assert p.doi, f"{p.name}: DOI yok"
        assert p.provenance, f"{p.name}: dizinin nereden geldiği yazılmamış"
        assert p.readout, f"{p.name}: hangi deney ölçtü yazılmamış"
        assert p.host_serum, f"{p.name}: konak yazılmamış"
        assert p.label in ("positive", "negative", "nonspecific", "mixed")


def test_set_negatif_de_iceriyor():
    """Yalnız pozitiflerden kurulu set özgüllük ölçemez."""
    assert len(lv.negatives()) >= 5
    assert len(lv.positives()) >= 5


def test_mixed_ve_nonspecific_pozitif_sayilmaz():
    for p in lv.all_peptides():
        if p.label in ("mixed", "nonspecific"):
            assert not p.is_positive, f"{p.name} pozitif sayılmamalı"


def test_diziler_tekil():
    seqs = [p.seq.upper() for p in lv.all_peptides()]
    assert len(seqs) == len(set(seqs)), "yinelenen dizi var"


def test_ozet_sayimlari_tutarli():
    s = lv.summary()
    assert s["n_total"] == len(lv.all_peptides())
    assert (s["n_positive"] + s["n_negative"] + s["n_nonspecific"]
            + s["n_mixed"]) == s["n_total"]
    assert s["n_studies"] >= 5
    assert s["scope_note"]                      # kapsam notu zorunlu
    assert s["missing_sequences"]               # dürüst eksik muhasebesi


def test_eslesme_siniflari():
    p = lv.positives()[0]
    hits = dict(lv.find_overlaps(p.seq))
    assert hits.get(p) == "exact"
    assert lv.find_overlaps("") == []


# --------------------------------------------------------------------------- #
# Korelat uyarıları
# --------------------------------------------------------------------------- #
def test_leptospira_korelat_uyarilari_zorunlu():
    assert correlate.is_mandatory("leptospira")
    ids = {c.id for c in correlate.caveats("Leptospira interrogans")}
    for need in ("lepto_lps_correlate_unmodelled",
                 "lepto_titre_does_not_predict",
                 "lepto_immunogenicity_can_harm"):
        assert need in ids, f"{need} eksik"


def test_modul_hicbir_organizmada_sessiz_kalmaz():
    cs = correlate.caveats("Bovine herpesvirus 1")
    assert cs, "tanınmayan organizmada da genel uyarı dönmeli"
    assert not correlate.is_mandatory("Bovine herpesvirus 1")
    assert correlate.caveats(None)


def test_critical_ve_major_uyarilar_atifli():
    for c in correlate.caveats("leptospira"):
        if c.severity in ("critical", "major"):
            assert c.citations, f"{c.id}: atıf yok"
            assert any("doi" in x.lower() for x in c.citations), f"{c.id}: DOI yok"


def test_uyarilar_iki_dilli_ve_siralanmis():
    cs = correlate.caveats("leptospira")
    sev = [c.severity for c in cs]
    assert sev == sorted(sev, key=lambda s: {"critical": 0, "major": 1}.get(s, 2))
    for c in cs:
        assert c.text("tr") and c.text("en") and c.text("tr") != c.text("en")
        assert c.title("tr") and c.title("en")


def test_profil_cozumleme():
    assert correlate.resolve("L. borgpetersenii") == "leptospira"
    assert correlate.resolve("leptospiroz aşısı") == "leptospira"
    assert correlate.resolve("Mycobacterium bovis") is None


# --------------------------------------------------------------------------- #
# Anchor: allele özgülük
# --------------------------------------------------------------------------- #
def test_mhc_ii_icin_varsayilan_pozisyon_uydurulmaz():
    """MHC-II çekirdek ofseti bilinmeden pozisyon verilmemeli."""
    assert mhc_motif.anchor_positions("mhc_ii", 15, allele=None) == []
    assert mhc_motif.peptide_anchors("A" * 15, "mhc_ii", allele=None) == {}


def test_mhc_i_arac_yoksa_dururst_varsayilan():
    """Allel bilinmiyorsa P2+PΩ'ye düşer (ve bu durum etiketlenir)."""
    assert mhc_motif.anchor_positions("mhc_i", 9, allele=None) == [2, 9]
    anch = mhc_motif.peptide_anchors("YLQPRTFLL", "mhc_i", allele=None)
    assert anch == {"P2": "L", "PΩ": "L"}


def test_varsayilan_yol_olculmedigini_soyler():
    class P:
        def __init__(s, seq):
            s.seq, s.metrics, s.methods = seq, {}, {}
    p = P("YLQPRTFLL")
    mhc_motif.annotate(p, "", "mhc_i", host="İnsan")
    # allel yok -> ölçüm yapılamaz; basis ya yoktur ya 'ÖLÇÜLMEDİ' der
    basis = p.metrics.get("anchor_basis", "")
    if basis:
        assert "ÖLÇÜLMEDİ" in basis
        assert "ALLELE ÖZGÜ" in basis.upper() or "allele özgü" in basis


def test_kisa_peptitte_anchor_yok():
    assert mhc_motif.anchor_positions("mhc_i", 2, allele=None) == []


def test_basis_satirlari_tekilleniyor():
    class P:
        def __init__(s, b):
            s.metrics = {"anchor_basis": b} if b else {}
    lines = mhc_motif.basis_lines([P("İnsan / HLA-A*02:01"), P("İnsan / HLA-A*02:01"),
                                   P("Sığır / BoLA-3:00101"), P(None)])
    assert lines == ["İnsan / HLA-A*02:01", "Sığır / BoLA-3:00101"]


def test_onbellek_surumu_eski_girdileri_gecersiz_kilar():
    """v1 (sabit P2+PΩ) girdileri v2'de kullanılmamalı."""
    assert mhc_motif._CACHE_VERSION >= 2
    assert mhc_motif._load_cache().get("_version") == mhc_motif._CACHE_VERSION


def test_esik_parametreleri_makul():
    assert 0 < mhc_motif.ANCHOR_SPREAD_FRAC < 1
    assert mhc_motif.ANCHOR_SPREAD_MIN > 0
    assert mhc_motif.MAX_ANCHORS >= 2


@pytest.mark.parametrize("cls,peplen,expect_last", [("mhc_i", 9, 9), ("mhc_i", 10, 10)])
def test_pOmega_peptit_uzunluguna_gore(cls, peplen, expect_last):
    pos = mhc_motif.anchor_positions(cls, peplen, allele=None)
    assert expect_last in pos, "C-terminal anchor peptidin sonuna bakmalı"


# --------------------------------------------------------------------------- #
# Rapor/arayüz bağlantısı
# --------------------------------------------------------------------------- #
def test_organizma_girdiden_cozulur():
    """Profil ('bacteria') organizmayı söylemez; girdi adından çözülmeli."""
    assert correlate.detect_from_input("data/uploads/leptospira_x.gb") == "leptospira"
    assert correlate.detect_from_input(None, extra="ORGANISM Leptospira borgpetersenii") == "leptospira"
    assert correlate.detect_from_input("data/uploads/bohv1.gb") is None


def test_rapor_korelat_bolumu_zorunlu(tmp_path):
    """HTML + PDF + JSON üçünde de uyarılar bulunmalı ve ÖZETTEN ÖNCE gelmeli."""
    import json as _json
    import re as _re

    from vaxforge import report
    from vaxforge.models import Peptide

    p = Peptide(seq="WIRVERMSAIMPDQIAKAAKAKPVQKL", parent="LipL32",
                kind="B-cell", start=134)
    p.passed = True
    meta = {"timestamp": "2026-11-23T00:00:00",
            "input": "leptospira_interrogans.gb",
            "input_path": "data/uploads/leptospira_interrogans.gb",
            "profile": "bacteria", "lang": "tr",
            "hosts": [{"name": "human", "label": "İnsan (Homo sapiens)",
                       "source": "IMGT/HLA", "n_mhc_i": 27, "n_mhc_ii": 15,
                       "predictors": {}}],
            "thresholds": [], "plan": [], "citations": []}
    paths = report.write_package(tmp_path, [p], meta)

    # JSON: makine-okunur
    j = _json.loads(Path(paths["json"]).read_text(encoding="utf-8"))
    corr = j["correlate"]
    assert corr["resolved"] == "leptospira"
    assert corr["mandatory"] is True
    assert corr["n_critical"] == 3
    assert corr["caveats"] and all(c["citations"] for c in corr["caveats"]
                                   if c["severity"] != "info")

    # HTML: bölüm var ve ÖZETTEN ÖNCE
    html = Path(paths["html"]).read_text(encoding="utf-8")
    heads = [m.group(1) for m in _re.finditer(r"<h2>(.*?)</h2>", html, _re.S)]
    assert heads, "rapor başlıksız"
    assert "KORELAT" in heads[0], f"korelat bölümü ilk değil: {heads[:2]}"
    assert "aglütine eden" in html          # anti-LPS korelatı
    assert "DAHA HIZLI" in html             # ağırlaşma bulgusu


def test_tanimayan_organizmada_da_bolum_yazilir(tmp_path):
    """Sessiz kalınmaz: genel uyarı her organizmada rapora girer."""
    from vaxforge import report
    from vaxforge.models import Peptide

    p = Peptide(seq="AAALAAAAV", parent="X", kind="B-cell", start=1)
    p.passed = True
    meta = {"timestamp": "2026-11-23T00:00:00", "input": "bohv1.gb",
            "input_path": "data/uploads/bohv1.gb", "profile": "virus",
            "lang": "tr", "hosts": [], "thresholds": [], "plan": [],
            "citations": []}
    paths = report.write_package(tmp_path, [p], meta)
    html = Path(paths["html"]).read_text(encoding="utf-8")
    assert "KORELAT" in html
    import json as _json
    corr = _json.loads(Path(paths["json"]).read_text(encoding="utf-8"))["correlate"]
    assert corr["resolved"] is None
    assert corr["n_caveats"] >= 1
    assert corr["mandatory"] is False
