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

# ══════════════════════════════════════════════════════════════════════
# 2026-11-12 bilimsel doğruluk düzeltmeleri — regresyon testleri
# Hepsi AĞSIZ ve HARİCİ ARAÇSIZ: yalnız config + saf Python mantığı.
# ══════════════════════════════════════════════════════════════════════

def test_localization_allowed_is_gram_aware():
    """Gram-pozitifte sitoplazmik zar DIŞ sınırdır → `membrane` izinli olmalı.

    Benchmark'ta iki koruyucu yüzey antijeni (CapB/B. anthracis, ActA/L.
    monocytogenes) tam bu yüzden elenmişti: PSORTb ikisine de `membrane` dedi,
    tek ortak bakteri listesi ise Gram-negatif mantığıyla yazılmıştı.
    Bkz. benchmark/MISS_ATTRIBUTION.md §2.
    """
    from vaxforge.config_loader import ThresholdConfig
    cfg = ThresholdConfig.load()

    pos = set(cfg.resolve("bacteria", gram="positive")["localization"].params["allowed"].value)
    neg = set(cfg.resolve("bacteria", gram="negative")["localization"].params["allowed"].value)
    unk = set(cfg.resolve("bacteria")["localization"].params["allowed"].value)

    # Gram-pozitif: zar yüzeydir; dış zar ve periplazma YOKTUR.
    assert "membrane" in pos
    assert "outer_membrane" not in pos and "periplasm" not in pos
    assert {"extracellular", "cell_wall"} <= pos

    # Gram-negatif: iç zar dış zarla perdeli → yüzeye açık değil, izinli DEĞİL.
    assert "membrane" not in neg
    assert {"extracellular", "outer_membrane", "cell_wall", "periplasm"} <= neg

    # Gram bilinmiyorsa eski ihtiyatlı liste korunur (geriye uyum).
    assert unk == neg

    # Sitoplazma hiçbir koşulda izinli değil.
    assert "cytoplasm" not in pos and "cytoplasm" not in neg


def test_recorded_psortb_membrane_antigens_pass_grampos_gate():
    """Kayıtlı GERÇEK PSORTb çıktısıyla kapı kararı: eskiden elenen, artık geçer.

    PSORTb 'Cytoplasmic Membrane' → `membrane` (bkz. psortb._LOC_MAP). Benchmark
    analizinde CapB ve ActA için kayıtlı çıktı: localization=membrane, TM=1.
    Bu test o kayıtlı etiketi kapıdan geçirir — PSORTb'yi yeniden koşmaz.
    """
    from vaxforge.config_loader import ThresholdConfig
    from vaxforge.psortb import _LOC_MAP
    cfg = ThresholdConfig.load()

    assert _LOC_MAP["cytoplasmicmembrane"] == "membrane"
    recorded_label = _LOC_MAP["cytoplasmicmembrane"]   # gerçek koşuda alınan etiket
    recorded_tm = 1                                     # gerçek koşuda alınan TM sayısı

    gp = cfg.resolve("bacteria", gram="positive")
    gn = cfg.resolve("bacteria", gram="negative")

    # Lokalizasyon kapısı: Gram-pozitifte GEÇER, Gram-negatifte ELENİR.
    assert recorded_label in set(gp["localization"].params["allowed"].value)
    assert recorded_label not in set(gn["localization"].params["allowed"].value)

    # TM kapısı ikisini de elemiyor (max_helices ≥ 1), yani tek engel lokalizasyondu.
    assert recorded_tm <= int(gp["transmembrane"].params["max_helices"].value)


def test_hard_filter_flags_are_actually_read():
    """`hard_filter` config'te ilan edilir ve funnel ONU okumalı.

    Önceden funnel lokalizasyon + konak-homolojisini koşulsuz sert uyguluyordu;
    YAML'da `hard_filter: false` yazmak hiçbir şeyi değiştirmiyordu. Varsayılanlar
    true, yani bu test varsayılan davranışın da korunduğunu sabitler.
    """
    from pathlib import Path

    from vaxforge.config_loader import ThresholdConfig
    cfg = ThresholdConfig.load()
    r = cfg.resolve("bacteria", gram="positive")

    # Varsayılanlar: ikisi de SERT (bilimsel duruş değişmedi).
    assert r["localization"].hard_filter is True
    assert r["human_homology"].hard_filter is True

    # Kod bayrağı gerçekten okuyor mu (sabit kodlanmış değil mi)? Kaynağı
    # DOSYADAN okuruz: funnel'ı import etmek ağır bağımlılık zinciri çeker,
    # bu test ise ağsız/araçsız kalmalı.
    src = (Path(__file__).resolve().parents[1] / "vaxforge" / "funnel.py").read_text()
    assert 'tools["localization"].hard_filter' in src
    assert 'tools["human_homology"].hard_filter' in src


def test_cai_excludes_single_codon_families_and_has_baseline():
    """Sharp & Li 1987: Met/Trp aileleri CAI'den hariç; taban CAI anlamlı olmalı."""
    from vaxforge import cloning

    # Tek-kodonlu aileler indeksten çıkarılmış.
    assert cloning._SINGLE_CODON == {"ATG", "TGG"}
    assert "ATG" not in cloning._W_CAI and "TGG" not in cloning._W_CAI

    # En-sık kodonla yazılan dizi → CAI tanım gereği 1.0 (totoloji).
    prot = "MKLVGGSAAYDEK"
    assert cloning.cai(cloning.back_translate(prot)) == pytest.approx(1.0, abs=1e-6)

    # Taban çizgisi ANLAMLI olmalı: 0 < taban < 1, yani optimizasyon kazancı ölçülebilir.
    base = cloning.expected_cai(prot)
    assert 0.0 < base < 1.0

    # Yalnız Met/Trp → CAI tanımsız (0.0 sentinel), "kötü adaptasyon" değil.
    assert cloning.cai(cloning.back_translate("MWMW")) == 0.0


def test_overlapping_n_sequons_are_counted():
    """Bitişik sequon'lar örtüşür; örtüşmeyen tarama onları eksik sayıyordu."""
    from vaxforge import glyco

    # "NNSS": hem 1 (NNS) hem 2 (NSS) sequon → 2 bulunmalı (önceden 1 bulunuyordu).
    sc = glyco.scan_sequons("NNSS")
    assert sc["n_count"] == 2
    assert [pos for pos, _ in sc["n_sequons"]] == [1, 2]

    # Pro kuralı korunuyor: N-P-x asla sequon değil.
    motifs = [m for _, m in glyco.scan_sequons("NGSXNPSNAT")["n_sequons"]]
    assert motifs == ["NGS", "NAT"]

def test_every_i18n_key_has_both_tr_and_en():
    """Proje kuralı: arayüz/rapor metinleri DAİMA TR + EN.

    Tek dilli bir anahtar eklenirse bu test düşer. Uyarı metinleri de bu kurala
    tabidir: tek dilli uyarı, EN raporda sessizce TR görünmesi demektir.
    """
    from vaxforge.i18n import STRINGS, t

    missing_en = sorted(k for k, v in STRINGS.items() if not (v.get("en") or "").strip())
    missing_tr = sorted(k for k, v in STRINGS.items() if not (v.get("tr") or "").strip())
    assert missing_en == [], f"EN karşılığı eksik: {missing_en}"
    assert missing_tr == [], f"TR karşılığı eksik: {missing_tr}"

    # 2026-11-12 düzeltmeleriyle gelen uyarıların ikisi de gerçekten farklı dilde.
    for key in ("cai_baseline", "cai_tautology_warn", "soft_localization",
                "soft_self_homology", "gram_allowed_note", "diag_unexpected"):
        assert key in STRINGS, key
        assert t("tr", key) != t("en", key), f"{key}: TR ve EN aynı metin"


def test_new_warnings_are_emitted_in_both_languages():
    """Yeni uyarılar veri katmanında da iki dilli taşınmalı (reason_en / note_en)."""
    from pathlib import Path as _P

    repo = _P(__file__).resolve().parents[1]

    # Huni yumuşak-mod izleri hem reason hem reason_en taşıyor.
    fsrc = (repo / "vaxforge" / "funnel.py").read_text()
    assert fsrc.count('"reason_en"') >= 2
    assert 'soft_localization' in fsrc and 'soft_self_homology' in fsrc

    # Klonlama notu iki dilli.
    from vaxforge import cloning
    csrc = (repo / "vaxforge" / "cloning.py").read_text()
    assert '"note_en"' in csrc
    assert callable(cloning.expected_cai)

    # Rapor klonlama bölümü yeni alanları ve uyarıyı i18n ile gösteriyor.
    rsrc = (repo / "vaxforge" / "report.py").read_text()
    for key in ("cai_optimized", "cai_baseline", "cai_tautology_warn"):
        assert key in rsrc, key

    # Tanılama bölümü de i18n ile basılıyor, yani EN raporda İngilizce görünür.
    assert "_diagnostics_section" in rsrc
    for key in ("diagnostics", "diag_clean", "diag_unexpected"):
        assert f'"{key}"' in rsrc or f"'{key}'" in rsrc, key

    # Eklenen her i18n anahtarı GERÇEKTEN kullanılıyor mu (ölü anahtar yok)?
    from vaxforge.i18n import STRINGS
    used = rsrc + fsrc + csrc + (repo / "app.py").read_text()
    new_keys = ["cai_optimized", "cai_baseline", "cai_skipped", "cai_tautology_warn",
                "cai_skipped_warn", "soft_localization", "soft_self_homology",
                "diagnostics", "diag_clean", "diag_unexpected", "gram_allowed_note"]
    unused = [k for k in new_keys if k in STRINGS and k not in used]
    assert unused == [], f"kullanılmayan i18n anahtarı: {unused}"

def test_immforge_availability_checks_entrypoint_not_just_directory():
    """`available()` dizin varlığına değil, GERÇEK giriş noktasına bakmalı.

    Eski davranış: ImmForge dizini varsa "araç var" denip adım çalıştırılıyor,
    `immforge.integrations.vaxforge` bulunamadığı için ModuleNotFoundError ile
    düşüyordu. Yani dürüst ARAÇ_YOK yerine BEKLENMEYEN istisna üretiliyordu.
    Bu test, yokluğun etiketli ve gerekçeli olmasını sabitler.
    """
    from vaxforge import immforge_score as im

    ok, why = im.entrypoint_available()
    assert isinstance(ok, bool) and isinstance(why, str) and why

    # Giriş noktası yoksa available() False olmalı (dizin var olsa bile).
    if not ok:
        assert im.available() is False
        from pathlib import Path as _P
        assert _P(im.IMMFORGE_PATH).is_dir() or "dizini yok" in why

        class _Host:
            name = "human"

        r = im.score_construct("MKLVGGSAAYDEK", [_Host()])
        assert r["verdict"] == "ARAÇ_YOK"
        assert r.get("reason")          # gerekçe taşınıyor
        assert r["tool"] == "none"

    # Kalibrasyon katmanı yol çözümüyle bulunuyor (ayrı ve bağımsız yetenek).
    assert isinstance(im.calibrated_species(), set)

def test_unconfigured_dtu_wrapper_is_not_reported_as_runnable(tmp_path):
    """Yapılandırılmamış DTU sarmalayıcısı 'çalışır' sayılmamalı.

    DTU paketleri `setenv NMHOME /tools/src/...` gibi SATICI YER TUTUCU yoluyla
    gelir; kurulumda elle düzeltilmezse araç hiç hesaplama yapmaz. Eski kontrol
    yalnız 'format error' metnine bakıyordu, bu yüzden `runnable()` True dönüyor
    ve `predict()` sessizce heuristik vekile düşüyordu. Sonuç: vekil sayılarla
    tamamlanan bir koşu GERÇEK araç koşusu sanılabiliyordu.
    """
    from vaxforge.netmhc_local import _wrapper_home_ok

    # Yer tutucu yol: dizin yok → çalışmaz sayılmalı.
    bad = tmp_path / "netMHCpan"
    bad.write_text("#!/bin/tcsh\nsetenv\tNMHOME\t/tools/src/netMHCpan-4.2\n")
    assert _wrapper_home_ok(str(bad)) is False

    # Gerçek dizini gösteren yol → çalışır sayılmalı.
    home = tmp_path / "netMHCpan-4.2"
    home.mkdir()
    good = tmp_path / "netMHCpan_ok"
    good.write_text(f"#!/bin/tcsh\nsetenv\tNMHOME\t{home}\n")
    assert _wrapper_home_ok(str(good)) is True

    # Bu deseni hiç kullanmayan araç (TMHMM/SignalP) hakkında karar VERMEMELİ.
    neutral = tmp_path / "tmhmm"
    neutral.write_text("#!/usr/bin/perl\nprint 1;\n")
    assert _wrapper_home_ok(str(neutral)) is True

    assert _wrapper_home_ok(None) is False



def test_mhc_i_allele_normalization_matches_tool_enforced_format():
    """SLA/BoLA tipi adlarda '*' -> ':' olmalı, HLA tipinde '*' atılmalı.

    NEDEN: netMHCpan '-a' doğrulamasını data/allelenames dosyasına karşı yapar.
    'SLA-1*0101' adından '*' atınca 'SLA-10101' çıkar; bu ad '-listMHC'
    çıktısında eş-anlamlı olarak GÖRÜNÜR ama allelenames'te YOKTUR → çağrı
    komple düşer, konak heuristik proxy'ye iner (kuduz koşusunda domuz 19 → 0
    bağlanan). Doğru biçim 'SLA-1:0101'.
    """
    from vaxforge.netmhc_local import _norm_allele
    # HLA: '*' sonrası kolon var → yalnız '*' atılır
    assert _norm_allele("HLA-A*01:01", "mhc_i") == "HLA-A01:01"
    assert _norm_allele("HLA-B*07:02", "mhc_i") == "HLA-B07:02"
    # SLA: '*' sonrası kolon YOK → '*' kolona çevrilir (haneler kaynaşmaz)
    assert _norm_allele("SLA-1*0101", "mhc_i") == "SLA-1:0101"
    assert _norm_allele("SLA-3*0401", "mhc_i") == "SLA-3:0401"
    # Zaten kolonlu / yıldızsız adlar dokunulmadan geçer
    assert _norm_allele("BoLA-1:00901", "mhc_i") == "BoLA-1:00901"
    assert _norm_allele("H-2-Kb", "mhc_i") == "H-2-Kb"


def test_allele_validation_uses_enforced_list_not_cli_listing():
    """Doğrulama '-listMHC' ÜST kümesine değil, zorlayıcı veri dosyasına yapılmalı.

    '-listMHC' eş-anlamlı yazımları da basar; '-a' bayrağı yalnız
    data/allelenames (MHC-I) ya da data/allelelist.txt (MHC-II) içeriğini kabul
    eder. Geniş listeye doğrulama yapmak geçersiz alleli süzgeçten geçirir.
    """
    import inspect

    from vaxforge import netmhc_local
    src = inspect.getsource(netmhc_local._supported)
    # Zorlayıcı liste CLI listesinden ÖNCE denenmeli. Docstring de '-listMHC'
    # içerdiği için konum karşılaştırması KOD satırlarına yapılır.
    i_enf = next(i for i, ln in enumerate(src.splitlines())
                 if "_enforced_allele_names(" in ln and "def " not in ln)
    i_cli = next(i for i, ln in enumerate(src.splitlines())
                 if "flag =" in ln)
    assert i_enf < i_cli, "zorlayıcı liste CLI listesinden sonra okunuyor"
    enf = inspect.getsource(netmhc_local._enforced_allele_names)
    assert "allelenames" in enf and "allelelist.txt" in enf


def test_configured_host_alleles_all_pass_enforced_validation():
    """config/hosts.yaml'daki her allel, kurulu araçların kabul ettiği biçimde olmalı.

    Araçlar kurulu değilse test atlanır (CI'da hafif kalsın).
    """
    import pytest
    import yaml

    from vaxforge import netmhc_local
    sup_i = netmhc_local._supported("mhc_i")
    sup_ii = netmhc_local._supported("mhc_ii")
    if not sup_i and not sup_ii:
        pytest.skip("netMHCpan/netMHCIIpan bu ortamda kurulu değil")
    cfg = yaml.safe_load(open("config/hosts.yaml"))
    hosts = cfg.get("hosts", cfg)
    bad = []
    for hname, h in hosts.items():
        if not isinstance(h, dict):
            continue
        for cls, sup in (("mhc_i", sup_i), ("mhc_ii", sup_ii)):
            if not sup:
                continue
            for a in (h.get(cls) or []):
                if netmhc_local._norm_allele(a, cls) not in sup:
                    bad.append(f"{hname}/{cls}/{a}")
    assert not bad, f"araçça tanınmayan allel(ler): {bad}"


def test_citations_are_marked_when_tool_is_unavailable():
    """Atıf tablosu kurulu OLMAYAN aracı 'kullanıldı' gibi göstermemeli.

    NEDEN: `citations.for_report()` eskiden KOŞULSUZ dönüyordu. NetCTL binary'si
    segfault verdiği koşuda 85 adayın tamamında `processing_norm` None kaldı ama
    rapor NetCTL'i Larsen 2005 atfıyla yöntem olarak listeledi.
    """
    from vaxforge import citations
    fake = {
        "netctl": {"label": "NetCTL-1.2", "available": False},
        "netmhc_local": {"label": "NetMHCpan / NetMHCIIpan", "available": True},
    }
    rows = citations.for_report(fake)
    by_tool = {r["tool"]: r for r in rows}
    assert by_tool["NetCTL-1.2 (kesim + TAP)"]["used"] is False
    assert by_tool["NetMHCpan-4.1/4.2"]["used"] is True
    # eşlenmemiş satır iddia ETMEZ (None), False değil
    unmapped = [r for r in rows if r["tool"] not in citations._TOOL_TO_STATUS]
    assert unmapped and all(r["used"] is None for r in unmapped)
    # durum verilmezse hiçbir satır 'kullanıldı' iddia etmez
    assert all(r["used"] is None for r in citations.for_report())


def test_tool_status_records_every_probe_without_raising():
    """tool_status.collect() hiçbir aracın available() istisnasını YUTMAMALI."""
    from vaxforge import tool_status
    st = tool_status.collect()
    assert len(st) == len(tool_status._PROBES)
    for key, entry in st.items():
        assert isinstance(entry["available"], bool), key
        assert entry["label"], key
        # istisna olduysa dürüstçe not düşülmüş ve available False olmalı
        if entry["note"]:
            assert entry["available"] is False, key
    s = tool_status.summarize(st)
    assert s["available"] + s["missing"] == s["total"] == len(st)


def test_netctl_availability_is_smoke_tested_not_just_file_existence():
    """netCTL sarmalayıcı dosyası varken bile binary çöküyorsa available() False.

    netCTL-1.2a'nın x86_64 binary'si bu çekirdekte (satıcının kendi test.fsa
    dosyasıyla bile) segfault veriyor.
    """
    import inspect

    from vaxforge import netctl
    src = inspect.getsource(netctl.available)
    assert "subprocess.run" in src, "available() gerçek bir duman testi yapmıyor"
    assert "returncode" in src


def test_bepipred_availability_validates_wrapper_home():
    """BepiPred sarmalayıcısındaki `setenv BP` yolu doğrulanmalı.

    Yer tutucu yol kalırsa araç `cannot find BP in ...` diyip çıkar,
    predict_residues BOŞ liste döner ve çağıran bunu 'B-hücre epitopu yok'
    sanar — ayrıca çalışan Kolaskar-Tongaonkar yedeğini de DEVRE DIŞI bırakır.
    """
    import inspect

    from vaxforge import bepipred
    src = inspect.getsource(bepipred._configured)
    assert "setenv" in src and "BP" in src
    assert "is_dir()" in src


def test_population_coverage_does_not_claim_real_when_tool_absent():
    """Araç yoksa yöntem etiketi 'GERÇEK (...)' dememeli."""
    from vaxforge import population
    assert population.METHOD.startswith("GERÇEK")
    assert population.METHOD_UNAVAILABLE.startswith("ARAÇ_YOK")
    import inspect
    src = inspect.getsource(population)
    assert "METHOD if _ok else METHOD_UNAVAILABLE" in src


def test_tool_status_labels_are_bilingual():
    """Araç etiketleri EN raporda Türkçe görünmemeli (proje kuralı: daima TR+EN)."""
    from vaxforge import tool_status
    assert all(len(p) == 3 for p in tool_status._PROBES), "etiketler TR+EN olmalı"
    st = tool_status.collect()
    for key, e in st.items():
        assert e.get("label") and e.get("label_en"), key
    summ = tool_status.summarize(st)
    assert "missing_labels" in summ and "missing_labels_en" in summ
    assert len(summ["missing_labels"]) == len(summ["missing_labels_en"])
    # Türkçeye özgü harf EN listesinde olmamalı
    tr_only = set("çğıöşüÇĞİÖŞÜ")
    bad = [x for x in summ["missing_labels_en"] if tr_only & set(x)]
    assert not bad, f"EN etiketinde Türkçe karakter: {bad}"


def test_toxicity_is_not_hard_filtered_outside_validated_length_domain():
    """ToxinPred2 8-15 aa epitoplarda SERT FİLTRE olarak kullanılmamalı.

    Kanıt: aracın altı referans setinin tamamında en kısa dizi 35 aa. Bu
    ortamda ölçüldü — tam dizide AUC 1.00, 9-mer'de 0.69; varsayılan eşik 0.6'da
    toksik-OLMAYAN proteinden gelen 9-mer'lerin 66/100'ü "toksik" çıkıyor. Araç
    bozukken bu kalibrasyon sorunu görünmüyordu; düzeltilince güvenli
    epitopların üçte ikisi sessizce elenecekti.
    """
    import inspect

    from vaxforge import survival, toxinpred
    assert toxinpred.MIN_VALIDATED_LENGTH == 35
    assert not toxinpred.in_validated_domain("SEFREAILL")      # 9 aa
    assert toxinpred.in_validated_domain("A" * 40)
    src = inspect.getsource(survival.run) if hasattr(survival, "run") \
        else inspect.getsource(survival)
    assert "in_validated_domain" in src, "huni uzunluk alanını kontrol etmiyor"
    # alan dışı sayacı özete yazılmalı
    assert "toksisite_alan_disi" in src


def test_out_of_domain_toxicity_warning_is_bilingual():
    """Alan dışı uyarısı TR + EN olmalı (proje kuralı)."""
    from vaxforge.i18n import STRINGS, t
    for key in ("tox_out_of_domain", "tox_ood_count"):
        assert STRINGS[key]["tr"].strip() and STRINGS[key]["en"].strip(), key
    assert "35" in t("tr", "tox_out_of_domain")
    assert "35" in t("en", "tox_out_of_domain")


def test_tool_status_uses_deep_verify_when_module_offers_it():
    """`available()` yalan söyleyebilir → `verify()` varsa GERÇEK koşu denenmeli.

    TMHMM (bozuk shebang), IApred (eksik bağımlılık) ve ToxinPred2 (eski
    sklearn pickle + pandas 1.x API) üçü de `available()` True derken hiç
    çalışmıyordu; huni onları "GERÇEK" olarak raporluyordu.
    """
    from vaxforge import iapred, tmhmm_local, tool_status, toxinpred
    for mod in (tmhmm_local, iapred, toxinpred):
        assert hasattr(mod, "verify"), mod.__name__
    st = tool_status.collect()
    for key in ("tmhmm_local", "iapred", "toxinpred"):
        assert st[key]["check"] == "verified", key
    s = tool_status.summarize(st)
    assert s["verified"] + s["shallow"] == s["total"]
    assert s["verified"] >= 3


def test_tool_check_depth_note_is_bilingual_and_formats():
    """Kontrol derinliği notu TR + EN olmalı ve sayıları yerleştirmeli."""
    from vaxforge.i18n import t
    for lang in ("tr", "en"):
        msg = t(lang, "tool_check_note").format(v=3, s=19)
        assert "3" in msg and "19" in msg
        assert "{v}" not in msg and "{s}" not in msg


def test_toxinpred_handles_single_peptide_input():
    """Satıcı script'i TEK dizi girdide düşüyordu (np.loadtxt 1B dizi).

    Araç kurulu değilse test atlanır.
    """
    import pytest

    from vaxforge import toxinpred
    if not toxinpred.available():
        pytest.skip("ToxinPred2 kurulu değil")
    assert len(toxinpred.predict(["SEFREAILLSEFRE"])) == 1


def test_psortb_gram_positive_membrane_antigens_pass_with_real_tool():
    """ActA ve CapB GERÇEK PSORTb'de `CytoplasmicMembrane` → Gram+ kapısından geçmeli.

    2026-11-13'te podman image'ı `--storage-opt ignore_chown_errors=true` ile
    çekilerek PSORTb ilk kez bu ortamda koşturuldu ve `MISS_ATTRIBUTION.md`
    kaydı BAĞIMSIZ olarak doğrulandı:

        P33379 ActA (L. monocytogenes) → CytoplasmicMembrane, skor 10.0
        P19580 CapB (B. anthracis)     → CytoplasmicMembrane, skor 10.0

    Gram-pozitifte sitoplazmik membran DIŞ sınırdır (dış membran yok), bu yüzden
    `membrane` meşru bir yüzey etiketidir. Gram-negatifte dış membranla
    perdelendiği için izinli DEĞİLDİR. Bu test araç olmadan da koşar: gerçek
    araç çıktısı sabit olarak gömülü.
    """
    from vaxforge.config_loader import ThresholdConfig
    from vaxforge.psortb import _LOC_MAP
    assert _LOC_MAP["cytoplasmicmembrane"] == "membrane"
    cfg = ThresholdConfig.load()
    pos = cfg.resolve("bacteria", gram="positive")["localization"].params["allowed"].value
    neg = cfg.resolve("bacteria", gram="negative")["localization"].params["allowed"].value
    none = cfg.resolve("bacteria", gram=None)["localization"].params["allowed"].value
    assert "membrane" in pos, "Gram+ membrane etiketi elenirse ActA/CapB kaybedilir"
    assert "membrane" not in neg, "Gram- sitoplazmik membran dış membranla perdeli"
    assert "membrane" not in none, "Gram bilinmiyorsa TEMKİNLİ (Gram-) liste"
    assert "outer_membrane" not in pos, "Gram+ bakteride dış membran YOK"


def test_psortb_offers_deep_verification():
    """`podman image exists` yetmez: image dururken de koşmayabilir."""
    import inspect

    from vaxforge import psortb
    assert hasattr(psortb, "verify")
    src = inspect.getsource(psortb.verify)
    assert "predict(" in src and "localization" in src


def test_bcell_method_label_reaches_the_exported_candidate_row():
    """B epitobunun yöntem sütunu BOŞ kalmamalı.

    Dışa aktarım yalnız `mhc_score` yöntemini okuyordu; B-hücre epitoplarında
    `method` NaN oluyordu, yani run.json/CSV okuyan biri skorun GERÇEK
    BepiPred'den mi klasik Kolaskar-Tongaonkar+Parker yedeğinden mi geldiğini
    göremiyordu. BepiPred sarmalayıcısı yer tutucu yolda kalıp sessizce yedeğe
    düşerken bu boşluk tam olarak onu gizledi.
    """
    import inspect

    from vaxforge import report
    src = inspect.getsource(report)
    i = src.index('"method": method_label(')
    snippet = src[i:i + 200]
    assert "bcell_score" in snippet, "yöntem sütunu B epitobu yöntemini okumuyor"
    # her iki anahtar da denenmeli
    assert "mhc_score" in snippet
