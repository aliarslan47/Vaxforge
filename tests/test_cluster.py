"""Örtüşen aday kümeleme testleri."""
from vaxforge import cluster
from vaxforge.models import Peptide


def _p(seq, parent="P1", kind="B", start=0, candidacy=0.5):
    return Peptide(seq=seq, parent=parent, kind=kind, start=start, candidacy=candidacy)


def test_kayan_pencere_komsulari_tek_kumede():
    """BoHV-1'deki asıl desen: 1 aa kaymış komşular tek epitop bölgesi."""
    peps = [
        _p("ASPSPPGP", start=0, candidacy=0.90),
        _p("SPSPPGPD", start=1, candidacy=0.87),
        _p("PSPPGPDG", start=2, candidacy=0.85),
    ]
    kept, s = cluster.cluster(peps, min_overlap=4)
    assert len(kept) == 1
    assert kept[0].seq == "ASPSPPGP"          # en yüksek puanlı temsilci
    assert kept[0].metrics["cluster_size"] == 3
    assert s["n_before"] == 3 and s["n_after"] == 1 and s["n_absorbed"] == 2


def test_temsilci_en_yuksek_puanli_olur_girdi_sirasi_ne_olursa():
    peps = [
        _p("ASPSPPGP", start=0, candidacy=0.10),
        _p("SPSPPGPD", start=1, candidacy=0.99),
    ]
    kept, _ = cluster.cluster(peps, min_overlap=4)
    assert len(kept) == 1 and kept[0].seq == "SPSPPGPD"


def test_esik_altinda_ortusme_ayri_kalir():
    """3 aa örtüşen peptitler min_overlap=4'te AYRI aday."""
    peps = [_p("AAAAAAAA", start=0), _p("AAAAABBB", start=5)]   # örtüşme = 3
    kept, s = cluster.cluster(peps, min_overlap=4)
    assert len(kept) == 2 and s["n_absorbed"] == 0


def test_farkli_protein_ayri_kalir():
    peps = [_p("ASPSPPGP", parent="P1", start=0), _p("ASPSPPGP", parent="P2", start=0)]
    kept, _ = cluster.cluster(peps, min_overlap=4)
    assert len(kept) == 2


def test_farkli_tur_ayri_kalir():
    """Aynı bölgeden çıkan B ve MHC-I epitopları farklı biyolojik iddialar."""
    peps = [_p("ASPSPPGP", kind="B", start=0), _p("ASPSPPGP", kind="MHC-I", start=0)]
    kept, _ = cluster.cluster(peps, min_overlap=4)
    assert len(kept) == 2


def test_uyeler_kaybolmaz_temsilcide_durur():
    peps = [
        _p("ASPSPPGP", start=0, candidacy=0.90),
        _p("SPSPPGPD", start=1, candidacy=0.87),
    ]
    kept, _ = cluster.cluster(peps, min_overlap=4)
    mem = kept[0].metrics["cluster_members"]
    assert [m["peptide"] for m in mem] == ["SPSPPGPD"]
    assert mem[0]["start"] == 1


def test_kapali_oldugunda_hicbir_sey_degismez():
    peps = [_p("ASPSPPGP", start=0), _p("SPSPPGPD", start=1)]
    kept, s = cluster.cluster(peps, min_overlap=4, enabled=False)
    assert len(kept) == 2 and s["enabled"] is False
    assert "cluster_size" not in peps[0].metrics


def test_puan_sirasi_korunur():
    peps = [
        _p("AAAAAAAA", parent="P1", start=0, candidacy=0.9),
        _p("BBBBBBBB", parent="P2", start=0, candidacy=0.8),
        _p("CCCCCCCC", parent="P3", start=0, candidacy=0.7),
    ]
    kept, _ = cluster.cluster(peps, min_overlap=4)
    assert [p.candidacy for p in kept] == [0.9, 0.8, 0.7]


def test_zincir_ortusme_tek_kumede_toplanir():
    """A-B örtüşür, B-C örtüşür ama A-C örtüşmez: temsilci A ise C ayrı kalır.

    Bu bilinçli bir seçim: temsilci-merkezli kümeleme, zincirin tamamını tek
    adaya indirip uzak bir bölgeyi gizlemesin.
    """
    peps = [
        _p("AAAAAAAA", start=0, candidacy=0.9),   # 0-8
        _p("AAAABBBB", start=4, candidacy=0.8),   # 4-12, A ile 4 örtüşür
        _p("BBBBCCCC", start=8, candidacy=0.7),   # 8-16, A ile 0 örtüşür
    ]
    kept, _ = cluster.cluster(peps, min_overlap=4)
    assert [p.seq for p in kept] == ["AAAAAAAA", "BBBBCCCC"]


def test_bos_girdi():
    kept, s = cluster.cluster([], min_overlap=4)
    assert kept == [] and s["n_after"] == 0
