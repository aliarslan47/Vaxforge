"""Düşük karmaşıklık dışlaması + %rank->skor ölçeği.

İkisi de 2026-11-22'de gerçek netMHCpan ölçümüyle bulunan kusurların testi:
  * `AAAAARAA` bağlanma eşiğini geçiyor (ölçülen %rank 1.131 <= 2.0) ama aşı
    adayı olamaz → konstrükte konmamalı, ★ almamalı, gerekçesi görünmeli.
  * `mhc_score = 1 - rank/100` eşik geçen peptitleri 0.98-1.00'e sıkıştırıyordu
    → güçlü/sınırda bağlayıcı adaylık puanında ayırt edilemiyordu.
"""

from vaxforge import complexity, mev, mhc_rank
from vaxforge.models import Peptide


# ── Karmaşıklık ölçütleri ──────────────────────────────────────────────────
def test_homopolimer_dusuk_karmasiklik():
    for pep in ("AAAAARAA", "AAAAAAAAA", "GGGGGGGGG", "AAARAAAAR",
                "AAPAAPAAPA", "RAAAAAAAR"):
        assert complexity.is_low_complexity(pep), pep


def test_gercek_guclu_baglayici_elenmez():
    """Ölçümde en güçlü bağlayıcılar arasında olan peptitler GEÇMELİ.

    `AGAPAYVAA` (%rank 0.008) Ala'dan zengin ama homopolimer değil — kural
    kompozisyonu değil TEKRARI eler.
    """
    for pep in ("TEHLRAVTA", "AGAPAYVAA", "AGATAPLAA", "RGVPVPPSL",
                "FKKTRLETY", "AAFDHFIRL"):
        assert not complexity.is_low_complexity(pep), pep


def test_metrikler_ve_gerekce():
    m = complexity.metrics("AAAAARAA")
    assert m["max_run"] == 5 and m["distinct"] == 2
    assert m["entropy"] < complexity.MIN_ENTROPY
    assert m["low_complexity"] is True
    assert m["reason"] and "homopolimer" in m["reason"]
    assert complexity.metrics("TEHLRAVTA")["reason"] is None


def test_annotate_passed_a_dokunmaz():
    """`passed` bağlanma kararıdır; karmaşıklık AYRI bir ölçüt."""
    p = Peptide(seq="AAAAARAA", parent="x", kind="MHC-I")
    p.passed = True
    complexity.annotate(p)
    assert p.passed is True                      # bağlanma kararı korunur
    assert p.metrics["low_complexity"] is True
    assert p.metrics["low_complexity_reason"]
    assert any("düşük karmaşıklık" in n for n in p.notes)
    assert "complexity" in p.methods


# ── %rank -> skor ölçeği ───────────────────────────────────────────────────
def test_rank_skoru_ayirt_ediyor():
    """Eski `1-rank/100` ile ÖLÇÜLEN sorun: yayılım 0.0196. Yeni: ~0.99."""
    guclu = mhc_rank.rank_to_score(0.003, "mhc_i", 0.5, 2.0)
    sinirda = mhc_rank.rank_to_score(1.131, "mhc_i", 0.5, 2.0)   # AAAAARAA
    assert guclu == 1.0
    assert sinirda < 0.25
    assert guclu - sinirda > 0.5          # eskiden 0.0113 idi


def test_rank_skoru_monoton_ve_sinirli():
    prev = 1.1
    for r in (0.01, 0.05, 0.1, 0.5, 1.0, 1.5, 1.99, 2.0, 10.0):
        s = mhc_rank.rank_to_score(r, "mhc_i", 0.5, 2.0)
        assert 0.0 <= s <= 1.0
        assert s <= prev                   # rank arttıkça skor azalır
        prev = s
    assert mhc_rank.rank_to_score(2.0, "mhc_i", 0.5, 2.0) == 0.0
    assert mhc_rank.rank_to_score(0.049, "mhc_i", 0.5, 2.0) == 1.0


def test_rank_skoru_sinif_esiklerine_demirli():
    """MHC-II eşikleri farklı (2/10): aynı rank iki sınıfta aynı skoru VERMEZ."""
    r = 1.0
    assert mhc_rank.rank_to_score(r, "mhc_i") < mhc_rank.rank_to_score(r, "mhc_ii")
    assert mhc_rank.rank_to_score(0.2, "mhc_ii") == 1.0     # 2.0/10
    assert mhc_rank.rank_to_score(10.0, "mhc_ii") == 0.0


def test_rank_skoru_bozuk_girdide_sifir():
    assert mhc_rank.rank_to_score("yok", "mhc_i") == 0.0
    assert mhc_rank.rank_to_score(-1, "mhc_i") == 0.0
    assert mhc_rank.rank_to_score(1.0, "mhc_i", 2.0, 0.5) == 0.0   # strong>weak


# ── MEV kaseti ─────────────────────────────────────────────────────────────
def _pep(seq, cand, kind="MHC-I"):
    p = Peptide(seq=seq, parent="prot1", kind=kind)
    p.passed = True
    p.candidacy = cand
    complexity.annotate(p)
    return p


def test_mev_dusuk_karmasikligi_disliyor_ve_iz_birakiyor():
    """Yüksek candidacy'li poly-Ala kasete GİRMEMELİ, ama iz BIRAKMALI."""
    peps = [_pep("AAAAARAA", 0.99),            # en yüksek skor, ama poly-Ala
            _pep("TEHLRAVTA", 0.80),
            _pep("FKKTRLETY", 0.75)]
    excluded: list[dict] = []
    picked = mev._pick(peps, "MHC-I", 6, excluded)
    seqs = [p.seq for p in picked]
    assert "AAAAARAA" not in seqs
    assert seqs == ["TEHLRAVTA", "FKKTRLETY"]
    assert len(excluded) == 1
    assert excluded[0]["seq"] == "AAAAARAA"
    assert "homopolimer" in excluded[0]["reason"]


def test_mev_construct_disladigini_rapora_tasiyor():
    peps = [_pep("AAAAARAA", 0.99), _pep("TEHLRAVTA", 0.80),
            _pep("GPGPVKNATL", 0.70, "MHC-II"), _pep("AAAAAAAAAAAAAAA", 0.95, "MHC-II")]
    c = mev.build_construct(peps, top_i=6, top_ii=4, top_b=4)
    assert "AAAAARAA" not in c.seq
    assert "TEHLRAVTA" in c.seq
    kinds = {e["kind"] for e in c.excluded_low_complexity}
    assert kinds == {"MHC-I", "MHC-II"}
    assert len(c.excluded_low_complexity) == 2


# ── Çıktı dosyaları ────────────────────────────────────────────────────────
def test_fasta_dislanani_aday_olarak_yazmaz(tmp_path):
    """top_peptides.fasta sentez için kullanılıyor: poly-Ala 'cand' olmamalı,
    ama dosyadan sessizce de kaybolmamalı (EXCLUDED etiketiyle sonda)."""
    from vaxforge import report

    peps = [_pep("AAAAARAA", 0.99), _pep("TEHLRAVTA", 0.80),
            _pep("FKKTRLETY", 0.75)]
    report.write_package(tmp_path, peps, {"lang": "tr", "mev": None})
    fa = (tmp_path / "top_peptides.fasta").read_text()
    cand = [ln for ln in fa.splitlines() if ln.startswith(">cand")]
    exc = [ln for ln in fa.splitlines() if ln.startswith(">excluded")]
    assert len(cand) == 2 and len(exc) == 1
    assert "AAAAARAA" in fa                      # bilgi kaybolmuyor
    assert all("AAAAARAA" not in c for c in cand)
    assert "EXCLUDED=low_complexity" in exc[0]
    assert "homopolimer" in exc[0]


def test_csv_karmasiklik_sutunlarini_tasiyor(tmp_path):
    import csv

    from vaxforge import report

    peps = [_pep("AAAAARAA", 0.99), _pep("TEHLRAVTA", 0.80)]
    report.write_package(tmp_path, peps, {"lang": "tr", "mev": None})
    rows = list(csv.DictReader((tmp_path / "candidates.csv").open()))
    by = {r["peptide"]: r for r in rows}
    assert by["AAAAARAA"]["low_complexity"] == "True"
    assert "homopolimer" in by["AAAAARAA"]["low_complexity_reason"]
    assert by["TEHLRAVTA"]["low_complexity"] == "False"
    assert float(by["TEHLRAVTA"]["seq_entropy"]) > 1.5
