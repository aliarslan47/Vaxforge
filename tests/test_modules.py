"""Birim testler — 11-faz eklentilerinin çekirdek (saf-fonksiyon) mantığı.

Hızlı + deterministik + AĞSIZ: sadece dış-araç/ağ gerektirmeyen mantık test edilir
(sequon taraması, kodon/CAI, BLOSUM mimikri, RSA matematiği, pLDDT parse, kalibrasyon
dosya taraması...). Ağ/ağır-araç yolları uçtan-uca smoke + manuel canlı testlerle
doğrulandı; burada regresyon güvencesi için mantık kilitlenir.
"""

from __future__ import annotations

import pytest


# ── Faz 2: DeepLoc lokalizasyon eşleme ──────────────────────────────
def test_deeploc_norm_mapping():
    from vaxforge import deeploc
    assert deeploc._norm("Extracellular") == "extracellular"
    assert deeploc._norm("Cell membrane") == "plasma_membrane"
    assert deeploc._norm("Cytoplasm") == "cytoplasm"
    assert deeploc._norm("Zznonsense") == "unknown"
    assert isinstance(deeploc.available(), bool)


# ── Faz 3: yapı (pLDDT parse + ölçek) / RSA / konformasyonel B ───────
def test_structure_plddt_rescale_and_sanitize():
    from vaxforge import structure
    assert structure._sanitize("mk-kx*L") == "MKKL"

    # PDB sabit-genişlik: B-factor kolon 61-66 (0-tabanlı [60:66]). Doğru hizala.
    def _atom(bf):
        ln = list(" " * 66)
        ln[0:6] = "ATOM  "
        ln[12:16] = " CA "
        ln[60:66] = f"{bf:6.2f}"
        return "".join(ln)
    pdb = _atom(0.40) + "\n" + _atom(0.60) + "\n"   # 0-1 ölçek → 0-100'e çevrilmeli
    vals = structure._ca_plddt(pdb)
    assert vals == [40.0, 60.0]
    assert structure.fold("", ) is None            # boş dizi
    assert structure.fold("A" * (structure.MAX_LEN + 1)) is None  # çok uzun


def test_rsa_maxasa_and_summary():
    from vaxforge import rsa
    assert len(rsa._MAX_ASA) == 20
    s = rsa.surface_summary([{"rsa": 0.5, "exposed": True},
                             {"rsa": 0.1, "exposed": False}])
    assert s["n_res"] == 2 and s["n_exposed"] == 1 and s["percent_surface"] == 50.0
    assert rsa.rsa_from_pdb("") is None


def test_bcell_conf_patches():
    from vaxforge import bcell_conf
    seq = "DDDDDDDD"      # hidrofilik (Asp) — yüzey B eğilimli
    rsa_list = [{"rsa": 0.9, "exposed": True} for _ in seq]
    out = bcell_conf.predict(seq, rsa_list)
    assert out and out["n_patches"] >= 1
    assert out["patches"][0]["length"] >= bcell_conf.MIN_PATCH
    assert 0.0 <= bcell_conf._norm_parker("D") <= 1.0
    assert bcell_conf.predict("AAA", []) is None


# ── Faz 4: glikozilasyon sequon ─────────────────────────────────────
def test_glyco_sequon_scan():
    from vaxforge import glyco
    sc = glyco.scan_sequons("NGSXNPSNAT")   # NGS ✓, NPS ✗ (Pro), NAT ✓
    motifs = [m for _, m in sc["n_sequons"]]
    assert "NGS" in motifs and "NAT" in motifs
    assert all(not m.startswith("NP") for m in motifs)  # N-P-x asla
    assert glyco.available() is True

    class P:
        def __init__(s, seq, kind): s.seq, s.kind, s.metrics, s.methods = seq, kind, {}, {}
    ps = [P("NASTKLPQR", "B"), P("AAAAAAA", "B")]
    summ = glyco.annotate_peptides(ps)
    assert summ["n_bcell"] == 2 and summ["n_bcell_masked"] == 1


# ── Faz 5: kodon optimizasyonu + CAI ────────────────────────────────
def test_cloning_backtranslate_cai_gc():
    from vaxforge import cloning
    dna = cloning.back_translate("MKLV")
    assert len(dna) % 3 == 0 and dna.endswith(cloning._STOP)
    assert cloning.cai(dna) == pytest.approx(1.0, abs=1e-6)  # hep en-sık kodon → CAI 1
    r = cloning.run("MKLVGGSAAYDEK")
    assert r and 0 <= r["gc_percent"] <= 100 and r["dna_len_bp"] % 3 == 0
    assert r["suggested_vector"] == "pET-28a(+)"


# ── Faz 6: VaxiJen çapraz-kontrol ───────────────────────────────────
def test_vaxijen_local_score():
    from vaxforge import vaxijen
    r = vaxijen.local_score("MKKLLFAIPLVVPFYSHSTPYEQ")
    assert set(r) >= {"score", "antigenic", "method"}
    assert 0.0 <= r["score"] <= 1.0 and isinstance(r["antigenic"], bool)


# ── Faz 7: BLOSUM Hamming-1 kısmi mimikri ───────────────────────────
def test_mimicry_near_hit_conservative_only():
    from vaxforge import mimicry
    kset = frozenset(["MKKLLFAIV"])
    # I→V konservatif (BLOSUM 3) → yakalanmalı
    hit = mimicry._near_hit("MKKLLFAII", kset)
    assert hit and hit[0] == "MKKLLFAIV" and hit[4] >= 1.0
    # D→W (BLOSUM -4) konservatif değil → None
    kset2 = frozenset(["WKKLLFAII"])
    assert mimicry._near_hit("DKKLLFAII", kset2) is None


# ── Faz 8: suş çekme taxon normalizasyonu ───────────────────────────
def test_strain_fetch_taxon_norm():
    from vaxforge import strain_fetch
    assert strain_fetch._norm_taxon("NCBITaxon:821") == "821"
    assert strain_fetch._norm_taxon("abc") is None
    bad = strain_fetch.fetch_strains("not-a-taxon")
    assert bad["ok"] is False and bad["paths"] == []


# ── Faz 9: konak-DB taxid çözümü ────────────────────────────────────
def test_host_db_taxid_resolution():
    from vaxforge import host_db
    assert host_db._resolve_taxid("rabbit", None) == 9986
    assert host_db._resolve_taxid("x", 12345) == 12345
    assert host_db._resolve_taxid("bilinmeyen", None) is None
    assert host_db.db_path("human").name == "human.dmnd"


# ── Faz 10: dosya-tabanlı ImmForge kalibrasyonu ─────────────────────
def test_immforge_calibration_filebased():
    from vaxforge import immforge_score as im
    cs = im.calibrated_species()
    assert isinstance(cs, set)
    # kalibre bir tür kendini, kalibre olmayan 'human'-proxy'yi kullanmalı
    assert im._immforge_species("chicken") == "human"
    if "mouse" in cs:
        assert im._immforge_species("mouse") == "mouse"


# ── HLA-hastalık anotasyonu (patojenite katmanı) ────────────────────
def test_hla_disease_annotation():
    from vaxforge import hla_disease as hd
    assert hd.available() is True
    # bilinen doğrulanmış allel — advers kategori + DOI taşımalı
    ann = hd.annotate("HLA-B*57:01")
    assert ann and ann[0]["category"] == "adverse"
    assert ann[0]["doi"] == "10.1056/NEJMoa0706135" and ann[0]["pmid"]
    # her kayıt PMID + DOI taşımalı (uydurma-yok güvencesi)
    for lst in hd._load().values():
        for e in lst:
            assert e.get("pmid") and e.get("doi") and e.get("category") in hd.CATEGORY_TONE
    # eşleşmeyen allel boş; panel özeti kategori sayar
    assert hd.annotate("HLA-A*99:99") == []
    s = hd.summarize(["HLA-B*57:01", "HLA-B*53:01", "HLA-A*99:99"])
    assert s["n_adverse"] == 1 and s["n_protective"] == 1 and len(s["entries"]) == 2


def test_localization_allowed_is_gram_aware():
    """Gram(+) bakteride PSORTb 'membrane' (sitoplazmik zar) meşru yüzey etiketidir.

    Gram-negatif mantığıyla yazılmış izinli liste CapB (B. anthracis) ve ActA
    (L. monocytogenes) gibi gerçek yüzey antijenlerini eliyordu (ikisi de membrane
    skor 10.0, TM=1). Bkz. benchmark/MISS_ATTRIBUTION.md §2.
    """
    from vaxforge.config_loader import ThresholdConfig
    from vaxforge import psortb

    cfg = ThresholdConfig.load()
    assert cfg.validate() == [], "config yapısal olarak temiz olmalı"

    params = cfg.resolve("bacteria")["localization"].params
    gram_neg = set(params["allowed"].value)
    gram_pos = set(params["allowed_gram_positive"].value)

    # PSORTb 'Cytoplasmic Membrane' -> 'membrane' eşlemesi korunuyor
    assert psortb._LOC_MAP["cytoplasmicmembrane"] == "membrane"

    # asıl regresyon: Gram+ listesi 'membrane'i kabul eder, Gram− listesi etmez
    assert "membrane" in gram_pos, "Gram(+) yüzey antijenleri elenmemeli"
    assert "membrane" not in gram_neg, "Gram(−)'de iç membran elenmeye devam etmeli"

    # Gram+'da anlamsız bölmeler listede olmamalı (dış zar / periplazma yok)
    assert "outer_membrane" not in gram_pos and "periplasm" not in gram_pos

    # sitoplazma hiçbir Gram tipinde kabul edilmez (güvenlik/özgüllük korunuyor)
    assert "cytoplasm" not in gram_neg and "cytoplasm" not in gram_pos

    # funnel gerçekten gram'a göre listeyi seçiyor mu (imza/akış kontrolü)
    import inspect
    from vaxforge import funnel
    src = inspect.getsource(funnel.run)
    assert "allowed_gram_positive" in src, "funnel.run gram-duyarlı listeyi seçmeli"
    assert "gram" in inspect.signature(funnel.run).parameters


def test_toxicity_gate_is_soft_and_length_appropriate():
    """Toksisite kapısı epitop düzeyinde YUMUŞAK ve peptit-düzeyi araç kullanıyor.

    REGRESYON (2026-10-10): ToxinPred2 (tam-protein modeli) 9-mer'lere sert filtre
    olarak uygulanıyordu; kalibrasyonda insan self 9-mer'lerinin %75'i ve
    doğrulanmış IEDB epitoplarının %63'ü "toksik" işaretlendi (influenza M1 58-66
    GILGFVFTL dahil) → ESAT-6/L7-L12/ActA/LLO "aday peptit üretmedi" durumuna
    düşüyordu. Bu test o mimarinin geri gelmesini engeller.
    """
    import inspect
    from vaxforge.config_loader import ThresholdConfig
    from vaxforge import survival

    cfg = ThresholdConfig.load()
    assert cfg.validate() == []

    tox = cfg.resolve("bacteria")["toxicity"]
    assert tox.hard_filter is False, "epitop düzeyinde toksisite sert filtre OLMAMALI"
    assert "hard_discard_score" in tox.params, "yüksek-güven eleme eşiği tanımlı olmalı"
    assert "toxicity_penalty" in tox.params, "yumuşak ceza katsayısı tanımlı olmalı"
    soft = float(tox.params["max_toxicity_score"].value)
    hard = float(tox.params["hard_discard_score"].value)
    assert soft < hard, "işaretleme eşiği eleme eşiğinin altında olmalı"

    src = inspect.getsource(survival.run)
    # peptit düzeyinde ToxinPred3; ToxinPred2'ye DÜŞMEK yasak (uzunluk sınıfı yanlış)
    assert "toxinpred3.predict" in src
    assert "toxinpred.predict" not in src
    # gerçek araç yoksa oyuncak sezgisel SERT filtre olarak kullanılmamalı
    assert "_tox_score" in inspect.getsource(survival), "vekil anotasyon olarak kalmalı"
    assert "ATLANDI" in src, "araç yoksa kapı atlanmalı (sessizce elemek yok)"

    # scoring yumuşak cezayı uyguluyor mu
    from vaxforge import scoring
    assert "toxicity_penalty_factor" in inspect.getsource(scoring.score)


def test_allergen_criterion_matches_sequence_length():
    """Alerjenite ölçütü uzunluk sınıfına göre seçiliyor.

    FAO/WHO 6-mer kuralı kısa peptitlerde özgül (9-mer'de yanlış-pozitif ~%2-3)
    ama TAM PROTEİNDE rastgele non-alerjen proteinlerin neredeyse tamamını
    işaretliyor. Bu yüzden MEV konstruktu (bir protein) için Codex/FAO-WHO'nun
    birincil ölçütü (80 aa pencerede ≥%35 kimlik) kullanılır.
    """
    import inspect
    from vaxforge import allergen, mev

    assert hasattr(allergen, "predict_protein")
    assert allergen._WINDOW == 80 and allergen._MIN_ID == 35.0

    src = inspect.getsource(mev)
    assert "allergen.protein_available()" in src
    assert "predict_protein" in src, "konstrukt protein ölçütüyle değerlendirilmeli"

    # ölçüt uygulanamıyorsa dürüst 'bilinmiyor' (None) döner, False DEĞİL
    out = allergen.predict_protein("")
    assert out["allergen"] is None
