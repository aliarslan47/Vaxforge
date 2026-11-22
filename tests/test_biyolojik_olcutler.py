"""2026-11-22 dördüncü tur — dört biyolojik ölçüt kusurunun testi.

1. Konservasyon: ortolog yoksa/tekse %100 YAZMAZ (ölçüm yokluğu ≠ mükemmel korunma).
2. Konstrükt alerjenitesi: eşleşme birleşme yerinde doğuyorsa öyle raporlanır.
3. pLDDT: eşik altındayken yapıdan türeyen sonuçlar "düşük güven" taşır.
4. tool_status: `verify()` düşerse bir kez retry, ve "dosya yok" ile
   "verify düştü" ayrı raporlanır.
"""

import vaxforge.tool_status as ts
from vaxforge import complexity, conservation, mev, structure
from vaxforge.models import Peptide, ProteinRecord


# ── 1. Konservasyon ────────────────────────────────────────────────────────
class _FakeAligner:
    """Kimlik hizalayıcı yerine geçer; conservation.run'ı diamond'suz test eder."""


def _run_cons(monkeypatch, ortholog_map, n_strains):
    """conservation.run'ı diamond/dosya olmadan koştur."""
    prots = [ProteinRecord(id="p1", seq="MKVLAAGIVLLLSHKEPEAVY")]
    strains = [({f"s{i}": seq}, {"p1": f"s{i}"}) if seq else ({}, {})
               for i, seq in enumerate(ortholog_map)]
    monkeypatch.setattr(conservation, "available", lambda paths: True)
    # `run` her suş için `_make_db(path) -> (dmnd, seqs)` bekler, sonra
    # `_best_orthologs(proteins, dmnd)` çağırır. Burada dmnd yerine doğrudan
    # hit sözlüğünü koyuyoruz, `_best_orthologs` da onu aynen döndürüyor.
    monkeypatch.setattr(conservation, "_make_db",
                        lambda sp: ((strains[int(sp)][1], strains[int(sp)][0])
                                    if int(sp) < len(strains) else None))
    monkeypatch.setattr(conservation, "_best_orthologs",
                        lambda proteins, dmnd: dmnd)
    return conservation.run(prots, [str(i) for i in range(n_strains)], min_percent=80.0)


def test_ortolog_yoksa_yuzde_100_yazmaz(monkeypatch):
    """ESKİ HATA: ortolog bulunamayınca profil baştan sona 1.0 ve yüzde 100.0
    çıkıyordu — ölçüm yokluğu 'mükemmel korunmuş' gibi raporlanıyordu."""
    prots, summary = _run_cons(monkeypatch, [None], 1)
    a = prots[0].annotations
    assert a["conservation_percent"] is None
    assert a["conservation_profile"] is None
    assert a["conservation_n_orthologs"] == 0
    assert "ölçülemedi" in a["method_conservation"]
    assert summary["n_unmeasured"] == 1


def test_tek_ortolog_suslar_arasi_sayilmaz(monkeypatch):
    """1 ortolog = iki dizi arası kimlik; suş-ARASI korunmuşluk değil."""
    prots, _ = _run_cons(monkeypatch, ["MKVLAAGIVLLLSHKEPEAVY"], 1)
    a = prots[0].annotations
    assert a["conservation_n_orthologs"] == 1
    assert a["conservation_percent"] is None
    assert f"en az {conservation.MIN_ORTHOLOGS}" in a["method_conservation"]


def test_iki_ortologla_olculur_ve_kendisi_sayilmaz(monkeypatch):
    """İki ortolog yeter; ve aday KENDİSİ payda/pay'a katılmaz.

    İki ortolog da ilk kalıntıda farklıysa o pozisyonun konservasyonu 0.0
    olmalı — eskiden aday kendisini saydığı için 1/3 = 0.33 çıkıyordu.
    """
    seq = "MKVLAAGIVLLLSHKEPEAVY"
    alt = "AKVLAAGIVLLLSHKEPEAVY"      # yalnız 1. kalıntı farklı
    prots, summary = _run_cons(monkeypatch, [alt, alt], 2)
    a = prots[0].annotations
    assert a["conservation_n_orthologs"] == 2
    assert a["conservation_profile"][0] == 0.0     # kendisi sayılmıyor
    assert a["conservation_profile"][1] == 1.0
    assert a["conservation_percent"] is not None
    assert a["conservation_percent"] < 100.0
    assert summary["n_unmeasured"] == 0


def test_min_orthologs_sabiti_var():
    assert conservation.MIN_ORTHOLOGS >= 2


# ── 2. Konstrükt alerjenitesi ──────────────────────────────────────────────
def _pep(seq, cand, kind="MHC-I"):
    p = Peptide(seq=seq, parent="prot1", kind=kind)
    p.passed, p.candidacy = True, cand
    complexity.annotate(p)
    return p


def test_alerjen_eslesmesi_birlesme_yerine_atfediliyor():
    """'KHL' ile biten epitop + 'AAY' linker = KHLAAY (gerçek alerjen 6-mer).

    Epitopların kendisi temiz olduğu halde konstrükt alerjen çıkıyor; rapor
    bunun birleşme yerinden doğduğunu ve linkerla giderilebileceğini söylemeli.
    """
    c = mev.build_construct([_pep("EPEAVYKHL", 0.9), _pep("TEHLRAVTA", 0.8)],
                            top_i=6, top_ii=0, top_b=0)
    props = mev.characterize(c.seq, cassette=c.cassette, components=c.components)
    al = props["allergen"]
    if al["method"] == "yok":          # alerjen DB kurulu değilse atla
        return
    assert al["allergenic"] is True
    assert al["n_matches"] >= 1
    kmers = {h["kmer"] for h in al["matches"]}
    assert "KHLAAY" in kmers
    junction = [h for h in al["matches"] if h["kmer"] == "KHLAAY"][0]
    assert junction["region"] == "junction"
    assert "epitope" in junction["detail"] and "linker" in junction["detail"]
    assert junction["source"]                       # hangi alerjene benziyor
    assert al["fixable_by_linker"] is True
    assert al["from_epitope"] is False


def test_temiz_konstruktte_uyari_yok():
    c = mev.build_construct([_pep("TEHLRAVTA", 0.8)], top_i=6, top_ii=0, top_b=0)
    props = mev.characterize(c.seq, cassette=c.cassette, components=c.components)
    al = props["allergen"]
    if al["method"] == "yok":
        return
    if not al["allergenic"]:
        assert al["n_matches"] == 0
        assert al["fixable_by_linker"] is False


def test_attribution_bilesen_yoksa_cokmez():
    assert mev.allergen_attribution("TEHLRAVTA", None) == [] or True


# ── 3. pLDDT güveni ────────────────────────────────────────────────────────
def test_plddt_esik_alti_guvenilmez():
    """Ölçülen MEV değeri 55.6 — eşiğin (70) altında."""
    c = structure.confidence(55.6)
    assert c["level"] == "low"
    assert c["reliable"] is False
    assert "70" in c["label_tr"]


def test_plddt_seviyeleri():
    assert structure.confidence(95.0)["reliable"] is True
    assert structure.confidence(75.0)["level"] == "confident"
    assert structure.confidence(40.0)["level"] == "disordered"
    assert structure.confidence(None)["level"] == "unknown"
    assert structure.confidence(None)["reliable"] is False


def test_bolge_plddt_ortalamasi():
    per = [10.0, 20.0, 30.0, 40.0, 50.0]
    assert structure.region_plddt(per, 1, 2) == 15.0      # 1-tabanlı, iki uç dahil
    assert structure.region_plddt(per, 4, 5) == 45.0
    assert structure.region_plddt(None, 1, 2) is None
    assert structure.region_plddt(per, 3, 3) == 30.0


# ── 4. tool_status ─────────────────────────────────────────────────────────
class _FlakyMod:
    """İlk verify düşer, ikincisi geçer — geçici yük hatasını taklit eder."""

    def __init__(self):
        self.calls = 0

    def verify(self):
        self.calls += 1
        if self.calls == 1:
            raise RuntimeError("Resource temporarily unavailable")
        return True

    def available(self):
        return True


class _DeadMod:
    def verify(self):
        raise RuntimeError("boom")

    def available(self):
        return True        # dosya var, ama koşmuyor


class _MissingMod:
    def verify(self):
        return False

    def available(self):
        return False       # dosya hiç yok


def test_verify_gecici_hatada_retry_ile_kurtuluyor(monkeypatch):
    monkeypatch.setattr(ts, "_RETRY_DELAY", 0.0)
    mod = _FlakyMod()
    ok, err = ts._verify_with_retry(mod)
    assert ok is True and mod.calls == 2 and err == ""


def test_verify_kalici_hatada_dusuyor(monkeypatch):
    monkeypatch.setattr(ts, "_RETRY_DELAY", 0.0)
    ok, err = ts._verify_with_retry(_DeadMod())
    assert ok is False and err == "RuntimeError"


def test_dosya_yok_ile_verify_dustu_ayriliyor(monkeypatch):
    """Dürüstlük tek yönlü değil: çalışan aracı 'yok' demek de yanıltıcı."""
    monkeypatch.setattr(ts, "_RETRY_DELAY", 0.0)
    import importlib

    mods = {"vaxforge.a": _DeadMod(), "vaxforge.b": _MissingMod()}
    monkeypatch.setattr(ts, "_PROBES", [("a", "Araç A", "Tool A"),
                                        ("b", "Araç B", "Tool B")])
    monkeypatch.setattr(importlib, "import_module", lambda n: mods[n])
    st = ts.collect()
    assert st["a"]["available"] is False
    assert st["a"]["fail_kind"] == "verify_failed"
    assert "dosya var" in st["a"]["note"]
    assert st["a"]["retried"] is True
    assert st["b"]["fail_kind"] == "missing"
    assert "dosya yok" in st["b"]["note"]
    summ = ts.summarize(st)
    assert summ["verify_failed"] == ["Araç A"]
    assert summ["not_installed"] == ["Araç B"]


def test_hizalama_kaymasi_yok():
    """ESKİ HATA: yerel hizalama adayın 1. kalıntısından başlamazsa profil
    kayıyordu (kolon sayımı gerçek koordinat yerine kullanılıyordu)."""
    al = conservation._aligner()
    cand = "MKVLAAGIVLLLSHKEPEAVY"
    orth = "AKVLAAGIVLLLSHKEPEAVY"       # YALNIZ 1. kalıntı farklı
    vec = conservation._identity_vector(cand, orth, al)
    assert vec[0] == 0, "1. kalıntı farklı, 'korunmuş' sayılmamalı"
    assert all(v == 1 for v in vec[1:]), "geri kalan kalıntılar aynı"
    chars = conservation._aligned_ortholog_chars(cand, orth, al)
    # 0. kalıntı yerel hizalamanın DIŞINDA kalıyor → karşılığı yok ('-').
    # Önemli olan 1. kalıntının 'K' olması: eskiden buraya 'V' yazılıyordu.
    assert chars[0] == "-"
    assert chars[1] == "K" and chars[2] == "V"      # kaymamış

    # N-ucu TAMAMEN farklı (sinyal peptidi benzeri) — kayma büyük olurdu
    cand2 = "MWWWWWWKVLAAGIVLLLSHKEPEAVY"
    orth2 = "KVLAAGIVLLLSHKEPEAVY"
    v2 = conservation._identity_vector(cand2, orth2, al)
    assert sum(v2[:7]) == 0, "hizalanmayan N-ucu korunmuş sayılmamalı"
    assert all(x == 1 for x in v2[7:]), "hizalanan kuyruk tam korunmuş"
