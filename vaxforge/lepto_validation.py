"""Leptospira'ya özgü DENEYSEL peptit doğrulama seti.

NEDEN BU MODÜL VAR
------------------
Vaxforge'un recall benchmark'ının paydası normalde IEDB'dir (bkz. iedb_match.py).
Leptospira'da bu çalışmaz. 2026-11-22'de IEDB IQ-API canlı sorgulandı:

    NCBITaxon:171   (Leptospira, cins)        ->  2 yapı kaydı
    NCBITaxon:173   (L. interrogans)          -> 14 yapı kaydı (11'i non-negatif)
    NCBITaxon:29507 (L. kirschneri)           ->  1
    172/174/28182/28183/28184/1049759         ->  0

Toplam ~16 kayıt. Üstelik bunların bir kısmı tek bir TANISAL çalışmadan
(fay gösterimi peptitleri) gelir ve GroEL kaydının kaynağı bir insan HSP60
otoantikor çalışmasıdır (miyaljik ensefalomiyelit) — Leptospira immünolojisi
değildir. Yani IEDB paydası pratikte ~10'dur ve recall ölçümü anlamsızdır.

Bu modül, hakemli literatürden DENEYSEL olarak test edilmiş Leptospira
peptitlerini tek bir yerde, KANIT SINIFI ve KAYNAĞIYLA toplar.

DÜRÜSTLÜK KURALLARI (bu repodaki desene uygun)
----------------------------------------------
  * Hiçbir dizi tahmin edilmedi. Her dizi ya makale metninden ya IEDB
    kaydından birebir alındı; `provenance` alanı hangisi olduğunu söyler.
  * NEGATİF sonuçlar da sete dâhil: deneyde tanınmayan ya da özgül olmayan
    peptitler `label="negative"` / `"nonspecific"` olarak tutulur. Bir
    doğrulama seti yalnız pozitiflerden kurulursa özgüllük ölçülemez.
  * Etiketler koruma DEĞİL, antikor/T-hücre TANIMASIDIR. `readout` alanı
    hangi okuma olduğunu açıkça yazar. Leptospira'da tanıma ile koruma
    arasındaki kopukluk belgelidir (LipL32 paradoksu; bkz. correlate.py).
  * Koordinatlar UniProt O34094 (LipL32, 272 aa) dizisine karşı DOĞRULANDI;
    `coords` alanındaki aralıklar hesaplanmıştır, makaleden kopyalanmamıştır.

KAYNAKLAR (PubMed; DOI'ler aşağıda her kayıtta)
-----------------------------------------------
  Lin X, ve ark. Clin Vaccine Immunol. 2010;17(5):778-83.
      doi:10.1128/CVI.00405-09                        PMID 20237196
  Pissawong T, ve ark. J Microbiol Immunol Infect. 2020;53(1):11-22.
      doi:10.1016/j.jmii.2017.12.006                  PMID 29500044
  Ramli SR, ve ark. PLoS Negl Trop Dis. 2019;13(1):e0007131.
      doi:10.1371/journal.pntd.0007131                PMID 30677033
  Prasad M, ve ark. Microb Pathog. 2020;149:104407.
      doi:10.1016/j.micpath.2020.104407               PMID 32758519
  Mariano IHM, ve ark. PLoS One. 2023;18(2):e0281344.
      doi:10.1371/journal.pone.0281344                PMID 36745643
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

METHOD = "GERÇEK (hakemli literatür, deneysel peptit kaydı)"
VERSION = 1
COMPILED = "2026-11-22"

# LipL32 referans dizisi, UniProt'un KENDİ kaydından okunur (elle kopyalanmaz).
# Dosya: refdata/leptospira/O34094_LipL32.fasta — UniProt O34094, 272 aa,
# L. interrogans sg Icterohaemorrhagiae sv Lai 56601.
# Diziyi kaynağından okumak, 'koordinat/dizi uydurma' hatasını yapısal olarak
# imkânsız kılar; verify_coords() bu dosyaya karşı denetler.
LIPL32_FASTA = (
    Path(__file__).resolve().parent / "refdata" / "leptospira" / "O34094_LipL32.fasta"
)


@lru_cache(maxsize=1)
def lipl32_seq() -> str:
    """UniProt O34094 dizisi. Dosya yoksa boş döner (dürüst 'doğrulanamadı')."""
    if not LIPL32_FASTA.exists():
        return ""
    return "".join(
        line.strip()
        for line in LIPL32_FASTA.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith(">")
    )


@dataclass(frozen=True)
class ValPeptide:
    """Deneysel olarak test edilmiş tek bir Leptospira peptidi."""

    seq: str
    name: str                 # makalede kullanılan ad (ör. "Pep7", "SIR16-D1")
    protein: str              # kaynak protein
    uniprot: str | None       # biliniyorsa UniProt aksesyonu
    label: str                # positive | negative | nonspecific | mixed
    readout: str              # hangi deney ölçtü
    host_serum: str           # antikorun/hücrenin geldiği konak
    pmid: str
    doi: str
    provenance: str           # dizinin nereden alındığı (makale metni / IEDB)
    coords: str | None = None  # UniProt dizisinde 1-tabanlı aralık (hesaplandı)
    iedb_id: int | None = None
    note: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)

    @property
    def is_positive(self) -> bool:
        """Yalnız net pozitifler. 'mixed' ve 'nonspecific' pozitif SAYILMAZ."""
        return self.label == "positive"


# --------------------------------------------------------------------------- #
# Lin 2010 — LipL32 üzerinde öngörülen B+T birleşik epitopları
# Makale LipL21 ve LipL32 için toplam 6 epitop öngörüyor; IEDB yalnız LipL32
# kaynaklı 4'ünü taşıyor. Koordinat eşleştirmesi, makalenin "seçilmiş" ilan
# ettiği LipL32(133-160) ve LipL32(221-247) ile hangi dizilerin örtüştüğünü
# kesinleştirdi (aşağıdaki coords değerleri hesaplandı).
# --------------------------------------------------------------------------- #
_LIN2010 = [
    ValPeptide(
        seq="WIRVERMSAIMPDQIAKAAKAKPVQKL",
        name="LipL32(134-160)", protein="LipL32", uniprot="O34094",
        label="positive",
        readout="fay gösterimi + Western blot + ELISA; immünize BALB/c CD4+ T "
                "lenfositleri tarafından tanındı, baskın Th1",
        host_serum="fare (BALB/c, H-2d) + tavşan antiserumu",
        pmid="20237196", doi="10.1128/CVI.00405-09",
        provenance="IEDB structure_id 133740 (dizi birebir)",
        coords="134-160", iedb_id=133740,
        note="Makalenin antijen olarak SEÇTİĞİ iki LipL32 epitobundan biri "
             "(makalede LipL32_133-160 olarak anılıyor).",
        tags=("B-cell", "T-cell", "Th1"),
    ),
    ValPeptide(
        seq="VKGSFVASVGLLFPPGIPGVSPLIHS",
        name="LipL32(222-247)", protein="LipL32", uniprot="O34094",
        label="positive",
        readout="fay gösterimi + Western blot + ELISA; immünize BALB/c CD4+ T "
                "lenfositleri tarafından tanındı, baskın Th1",
        host_serum="fare (BALB/c, H-2d) + tavşan antiserumu",
        pmid="20237196", doi="10.1128/CVI.00405-09",
        provenance="IEDB structure_id 133732 (dizi birebir)",
        coords="222-247", iedb_id=133732,
        note="Makalenin antijen olarak SEÇTİĞİ ikinci LipL32 epitobu "
             "(makalede LipL32_221-247).",
        tags=("B-cell", "T-cell", "Th1"),
    ),
    ValPeptide(
        seq="VKTLLPYGSVINYYGYVKPG",
        name="LipL32(47-66)", protein="LipL32", uniprot="O34094",
        label="negative",
        readout="öngörüldü ve fay gösteriminde ifade edildi, ancak Western "
                "blot/ELISA ile antijen olarak SEÇİLMEDİ",
        host_serum="fare (BALB/c) + tavşan antiserumu",
        pmid="20237196", doi="10.1128/CVI.00405-09",
        provenance="IEDB structure_id 133733 (dizi birebir)",
        coords="47-66", iedb_id=133733,
        note="Öngörü pozitifi / deney negatifi — özgüllük ölçümü için değerli.",
        tags=("B-cell", "predicted-not-selected"),
    ),
    ValPeptide(
        seq="KKLLVRGLYRISFTTYK",
        name="LipL32(202-218)", protein="LipL32", uniprot="O34094",
        label="negative",
        readout="öngörüldü ve fay gösteriminde ifade edildi, ancak Western "
                "blot/ELISA ile antijen olarak SEÇİLMEDİ",
        host_serum="fare (BALB/c) + tavşan antiserumu",
        pmid="20237196", doi="10.1128/CVI.00405-09",
        provenance="IEDB structure_id 133618 (dizi birebir)",
        coords="202-218", iedb_id=133618,
        note="Öngörü pozitifi / deney negatifi.",
        tags=("B-cell", "predicted-not-selected"),
    ),
]

# --------------------------------------------------------------------------- #
# Pissawong 2020 — LipL32 immünodominans haritası (truncation parçaları)
# --------------------------------------------------------------------------- #
_PISSAWONG2020 = [
    ValPeptide(
        seq="IAKAAKAKPVQKLDDDDDGDDTYKEERHNKYNSLTRI",
        name="LipL32(148-184)", protein="LipL32", uniprot="O34094",
        label="positive",
        readout="rekombinant truncation; leptospiroz hasta serumu ve tavşan "
                "anti-Leptospira serumuyla IgG-DOMİNANT",
        host_serum="insan (konvalesan, eşli serum) + tavşan",
        pmid="29500044", doi="10.1016/j.jmii.2017.12.006",
        provenance="IEDB structure_id 974358 (dizi birebir)",
        coords="148-184", iedb_id=974358,
        note="Patojenik ve ara Leptospira türleri arasında ~%70 kimlik; "
             "yüzeye yüksek erişilebilir bölge.",
        tags=("B-cell", "IgG-dominant"),
    ),
    ValPeptide(
        seq="KEERHNKYNSLTRIKIPNPPKSFDDLKNIDTKKLLVRGLYRISF",
        name="LipL32(171-214)", protein="LipL32", uniprot="O34094",
        label="positive",
        readout="rekombinant truncation; hasta serumuyla IgM-DOMİNANT",
        host_serum="insan (akut)",
        pmid="29500044", doi="10.1016/j.jmii.2017.12.006",
        provenance="IEDB structure_id 975034 (dizi birebir)",
        coords="171-214", iedb_id=975034,
        note="IgM-dominant bölge; akut faz tanımada güçlü, IgG'de daha zayıf.",
        tags=("B-cell", "IgM-dominant"),
    ),
]

# --------------------------------------------------------------------------- #
# Ramli 2019 — ORFeome fay gösterimi; 18 tekil klonun 5'i titrasyon ELISA ile
# doğrulandı, 2'si bireysel serumlarda anlamlı ayırt etti.
# Kalan 13 klonun DİZİSİ YAYINLANMADI (yalnız lokus + pozisyon verildi) →
# sete alınmadı; bkz. MISSING_SEQUENCES.
# --------------------------------------------------------------------------- #
_RAMLI2019 = [
    ValPeptide(
        seq="NPTAEDWEVDLVDTGALTMTGGRLLRLKDQLSKETFMVTY",
        name="SIR16-D1", protein="Glukoz-1-fosfat sitidiltransferaz (rfbF)",
        uniprot=None, label="positive",
        readout="sentetik biyotinli peptit ELISA, 16 akut hasta vs 16 sağlıklı; "
                "AUC 0,85 (0,70-1,0), duyarlılık 0,88, özgüllük 0,81",
        host_serum="insan (akut, IgM)",
        pmid="30677033", doi="10.1371/journal.pntd.0007131",
        provenance="makale Tablo 3 (dizi birebir); IEDB structure_id 873187",
        iedb_id=873187,
        note="En iyi tek peptit. Hücre-içi protein olmasına rağmen hümoral "
             "yanıt uyarıyor (yazarların 'moonlighting' yorumu).",
        tags=("B-cell", "IgM", "diagnostic"),
    ),
    ValPeptide(
        seq="EFSKTIVEKANQFWMMVRGEGAYSKPTRISQFSIQGLMREEDVLKTS",
        name="SIR16-H1", protein="Hipotetik protein (L. interrogans sv Naam)",
        uniprot=None, label="positive",
        readout="sentetik peptit ELISA; AUC 0,78 (0,59-0,96), duyarlılık 0,94, "
                "özgüllük 0,69",
        host_serum="insan (akut, IgM)",
        pmid="30677033", doi="10.1371/journal.pntd.0007131",
        provenance="makale Tablo 3 (dizi birebir); IEDB structure_id 871977",
        iedb_id=871977,
        note="İşlevi bilinmeyen proteinden; ters vaksinolojinin anotasyona "
             "bağımlı olmadığına örnek.",
        tags=("B-cell", "IgM", "diagnostic"),
    ),
    ValPeptide(
        seq="KSSISLGRGQAASVLAVPIWGRMYNRFYGGQNYPSFGED",
        name="SIR16-A1", protein="PBP1A ailesi penisilin-bağlayan protein",
        uniprot=None, label="mixed",
        readout="sentetik peptit ELISA; tek başına anlamlı ayırt etmedi, ancak "
                "rLipL32 + D1 + H1 ile birlikte sınıflandırmayı AUC 0,98'e çıkardı",
        host_serum="insan (akut, IgM)",
        pmid="30677033", doi="10.1371/journal.pntd.0007131",
        provenance="makale Tablo 3 (dizi birebir); IEDB structure_id 872713",
        iedb_id=872713,
        note="Tek başına negatif, kombinasyonda bilgi taşıyor → ikili etiket "
             "bu peptidi yanlış cezalandırır.",
        tags=("B-cell", "combination-only"),
    ),
    ValPeptide(
        seq="FGWNRDHFYLDGFFGSY",
        name="SIR16-C1", protein="Permeaz",
        uniprot=None, label="negative",
        readout="sentetik peptit ELISA; hasta ve kontrol serumları arasında "
                "anlamlı fark YOK",
        host_serum="insan (akut, IgM)",
        pmid="30677033", doi="10.1371/journal.pntd.0007131",
        provenance="makale Tablo 3 (dizi birebir); IEDB structure_id 872128",
        iedb_id=872128,
        tags=("B-cell",),
    ),
    ValPeptide(
        seq="EEPIYGLTKPVTTFDPVYT",
        name="SIR16-E6", protein="Sterol desatüraz",
        uniprot=None, label="negative",
        readout="sentetik peptit ELISA; anlamlı fark YOK",
        host_serum="insan (akut, IgM)",
        pmid="30677033", doi="10.1371/journal.pntd.0007131",
        provenance="makale Tablo 3 (dizi birebir); IEDB structure_id 871965",
        iedb_id=871965,
        tags=("B-cell",),
    ),
]

# --------------------------------------------------------------------------- #
# Mariano 2023 — LruC lipoproteininden 7 peptit; hamster (düşük-LPS bakterinle
# bağışıklanmış) ve insan hasta serumlarıyla Dot blot + ELISA.
# Pep1 ve Pep2'de sentezi kolaylaştırmak için BİRER amino asit değiştirilmiş
# (makalede bildirilmiş) → native dizi değil, note alanında işaretli.
# --------------------------------------------------------------------------- #
_MARIANO2023 = [
    ValPeptide(
        seq="WAILVPGA", name="Pep7", protein="LruC", uniprot=None,
        label="positive",
        readout="Dot blot (bağışıklanmış hamster serumu) EN YÜKSEK sinyal; "
                "insan hasta serumlarının TAMAMI tanıdı, kontrol tanımadı; "
                "ELISA ile doğrulandı",
        host_serum="hamster (düşük-LPS bakterin) + insan (hasta)",
        pmid="36745643", doi="10.1371/journal.pone.0281344",
        provenance="makale metni (dizi birebir); IEDB structure_id 2187865",
        iedb_id=2187865,
        note="Setin en güçlü pozitifi. NetMHC'de DRB3_0101'e güçlü bağlanma "
             "(13,70 nM, %0,08 rank). ELISA tanıması MAT titresiyle KORELE DEĞİL.",
        tags=("B-cell", "T-cell", "MHC-II", "ELISA-validated"),
    ),
    ValPeptide(
        seq="GSIPFTYNTVGQT", name="Pep1", protein="LruC", uniprot=None,
        label="positive",
        readout="Dot blot; yalnız bağışıklanmış hamster serumu tanıdı "
                "(kontrol tanımadı) → özgül",
        host_serum="hamster (düşük-LPS bakterin)",
        pmid="36745643", doi="10.1371/journal.pone.0281344",
        provenance="makale metni (dizi birebir)",
        note="DİKKAT: sentez için bir amino asit DEĞİŞTİRİLMİŞ — native LruC "
             "dizisi değil. MHC-II karşılığı bulunamadı (yalnız MHC-I).",
        tags=("B-cell", "sequence-modified"),
    ),
    ValPeptide(
        seq="VGTIPLNLVVTD", name="Pep2", protein="LruC", uniprot=None,
        label="mixed",
        readout="Dot blot hamsterda özgül POZİTİF; insan hasta serumlarında "
                "DÜŞÜK sinyal ve kontrolde de sinyal → insanda özgül değil",
        host_serum="hamster (pozitif) + insan (negatif)",
        pmid="36745643", doi="10.1371/journal.pone.0281344",
        provenance="makale metni (dizi birebir)",
        note="DİKKAT: sentez için bir amino asit değiştirilmiş. Konak-bağımlı "
             "sonuç: hamsterda pozitif, insanda değil. MHC karşılığı yok.",
        tags=("B-cell", "host-dependent", "sequence-modified"),
    ),
    ValPeptide(
        seq="YSSSFILIIKKG", name="Pep4", protein="LruC (ve diğerleri)",
        uniprot=None, label="nonspecific",
        readout="Dot blot; bağışıklanmış VE kontrol hamster serumlarının ikisi "
                "de tanıdı → özgül olmayan etkileşim",
        host_serum="hamster",
        pmid="36745643", doi="10.1371/journal.pone.0281344",
        provenance="makale metni (dizi birebir)",
        note="BLASTp'de birden çok proteinle eşleşiyor; HPLC'de safsızlık. "
             "Özgüllük testi için değerli negatif kontrol.",
        tags=("nonspecific", "multi-match"),
    ),
    ValPeptide(
        seq="AEGVSTVAYEDLYPSA", name="Pep3", protein="LruC", uniprot=None,
        label="negative",
        readout="Dot blot; bağışıklanmış ve kontrol serumları arasında anlamlı "
                "fark YOK (MHC bağlanma öngörüsü iyi olmasına rağmen)",
        host_serum="hamster",
        pmid="36745643", doi="10.1371/journal.pone.0281344",
        provenance="makale metni (dizi birebir); IEDB structure_id 2187593",
        iedb_id=2187593,
        note="İyi MHC öngörüsü / deney negatifi — öngörü-deney kopukluğunun "
             "doğrudan örneği.",
        tags=("predicted-not-recognised",),
    ),
    ValPeptide(
        seq="TKTVSSSD", name="Pep5", protein="LruC", uniprot=None,
        label="negative",
        readout="Dot blot; anlamlı etkileşim YOK",
        host_serum="hamster",
        pmid="36745643", doi="10.1371/journal.pone.0281344",
        provenance="makale metni (dizi birebir); IEDB structure_id 2187843",
        iedb_id=2187843,
        tags=("predicted-not-recognised",),
    ),
    ValPeptide(
        seq="LGSYPYDIFIKVK", name="Pep6", protein="LruC", uniprot=None,
        label="negative",
        readout="Dot blot; anlamlı etkileşim YOK",
        host_serum="hamster",
        pmid="36745643", doi="10.1371/journal.pone.0281344",
        provenance="makale metni (dizi birebir)",
        tags=("predicted-not-recognised",),
    ),
]

# --------------------------------------------------------------------------- #
# Prasad 2020 — LIC11334'ten tek B-hücre epitobu
# --------------------------------------------------------------------------- #
_PRASAD2020 = [
    ValPeptide(
        seq="EVENPDPKVAQE", name="KNSMP01", protein="LIC11334",
        uniprot=None, label="positive",
        readout="tavşan antiserumu ile Western blot, ELISA ve dot blot'ta "
                "antijenik",
        host_serum="tavşan (peptitle bağışıklanmış)",
        pmid="32758519", doi="10.1016/j.micpath.2020.104407",
        provenance="makale metni (dizi birebir); IEDB structure_id 1308568",
        iedb_id=1308568,
        note="Biyotinli sentez formu SGSGEVENPDPKVAQEC; buradaki dizi native "
             "epitop çekirdeğidir (linker/Cys hariç).",
        tags=("B-cell",),
    ),
]

PEPTIDES: list[ValPeptide] = (
    _LIN2010 + _PISSAWONG2020 + _RAMLI2019 + _MARIANO2023 + _PRASAD2020
)

# Dizisi yayınlanmadığı için sete ALINAMAYAN kayıtlar — dürüst muhasebe.
MISSING_SEQUENCES = {
    "30677033": (
        "Ramli 2019: 18 tekil fay gösterimi klonundan yalnız 5'inin dizisi "
        "yayınlandı (Tablo 3). Kalan 13 klon için yalnız NCBI lokusu ve "
        "pozisyon aralığı verilmiş (Tablo 2); dizi türetmek mümkün ama bazı "
        "aralıklar tutarsız (ör. '17046', '0-80') → türetilmedi."
    ),
    "24467929": (
        "Deveson Lucas 2014: 238 rekombinant protein sığır IFN-γ paneli "
        "PROTEİN düzeyinde, peptit düzeyinde değil. LipL32 için T-hücre "
        "reaktif bölge aa 20-200 olarak daraltılmış (peptit dizisi yok). "
        "Bu kaynak labels/ katmanına, bu modüle değil."
    ),
}

# Bu set hangi SORUYA cevap verir, hangisine VERMEZ — rapora yazılır.
SCOPE_NOTE = (
    "Bu doğrulama seti ANTİKOR/T-HÜCRE TANIMASINI ölçer, KORUMAYI ölçmez. "
    "Leptospira'da tanıma ile koruma arasındaki kopukluk belgelidir: LipL32 "
    "dış membran proteininin %75'ine kadarını oluşturur ve doğal enfeksiyonda "
    "baskın immünojendir, ancak lipL32 mutantı hamsterda akut hastalıkta ve "
    "sıçanda kronik kolonizasyonda yabanıl tipten ayırt edilemez "
    "(Murray 2009, doi:10.1128/IAI.01370-08)."
)


# --------------------------------------------------------------------------- #
# Erişim yardımcıları
# --------------------------------------------------------------------------- #
def all_peptides() -> list[ValPeptide]:
    return list(PEPTIDES)


def positives() -> list[ValPeptide]:
    """Net deneysel pozitifler (recall paydası)."""
    return [p for p in PEPTIDES if p.is_positive]


def negatives() -> list[ValPeptide]:
    """Deneyde tanınmayan ya da özgül olmayanlar (özgüllük/FP ölçümü)."""
    return [p for p in PEPTIDES if p.label in ("negative", "nonspecific")]


_HOST_KEYS = {
    "insan": "human", "human": "human",
    "fare": "mouse", "balb/c": "mouse", "c3h": "mouse",
    "hamster": "hamster",
    "tavşan": "rabbit", "rabbit": "rabbit",
    "sığır": "bovine", "bovine": "bovine",
}


def hosts_of(p: ValPeptide) -> list[str]:
    """Bu peptidin KANITININ geldiği konak(lar), Vaxforge konak adlarıyla.

    Bir peptit birden çok konakta test edilmiş olabilir (ör. Mariano 2023 Pep7:
    hamster + insan). Etiket konağa bağlıdır: aynı peptit bir konakta özgül
    pozitif, başka konakta özgül olmayabilir (Pep2 buna örnek).
    """
    s = (p.host_serum or "").lower()
    found = []
    for token, key in _HOST_KEYS.items():
        if token in s and key not in found:
            found.append(key)
    return found or ["unknown"]


def by_host() -> dict[str, list[ValPeptide]]:
    """Konak -> o konakta deneysel kanıtı olan peptitler."""
    out: dict[str, list[ValPeptide]] = {}
    for p in PEPTIDES:
        for h in hosts_of(p):
            out.setdefault(h, []).append(p)
    return out


def positives_for_host(host: str) -> list[ValPeptide]:
    """Belirli bir konak için recall paydası.

    Vaxforge çok-konaklı çalıştığında recall, KONAĞA GÖRE ölçülmelidir: insan
    serumuyla doğrulanmış bir epitobu fare konağı taramasında kaçırmak bir
    başarısızlık değildir. Payda raporda açıkça yazılır.
    """
    h = (host or "").lower()
    return [p for p in PEPTIDES if p.is_positive and h in hosts_of(p)]


def host_coverage() -> dict[str, dict[str, int]]:
    """Konak başına pozitif/negatif sayımı — raporun payda tablosu."""
    out: dict[str, dict[str, int]] = {}
    for host, peps in by_host().items():
        out[host] = {
            "n_total": len(peps),
            "n_positive": len([p for p in peps if p.is_positive]),
            "n_negative": len([p for p in peps
                               if p.label in ("negative", "nonspecific")]),
        }
    return out


def by_label() -> dict[str, list[ValPeptide]]:
    out: dict[str, list[ValPeptide]] = {}
    for p in PEPTIDES:
        out.setdefault(p.label, []).append(p)
    return out


def summary() -> dict:
    """Rapora ve testlere giden sayım özeti."""
    lab = by_label()
    return {
        "version": VERSION,
        "compiled": COMPILED,
        "method": METHOD,
        "n_total": len(PEPTIDES),
        "n_positive": len(lab.get("positive", [])),
        "n_negative": len(lab.get("negative", [])),
        "n_nonspecific": len(lab.get("nonspecific", [])),
        "n_mixed": len(lab.get("mixed", [])),
        "n_studies": len({p.pmid for p in PEPTIDES}),
        "pmids": sorted({p.pmid for p in PEPTIDES}),
        "n_with_iedb_id": len([p for p in PEPTIDES if p.iedb_id is not None]),
        "proteins": sorted({p.protein for p in PEPTIDES}),
        "hosts": sorted({p.host_serum for p in PEPTIDES}),
        # Recall KONAĞA GÖRE ölçülmeli: insan serumuyla doğrulanmış bir epitobu
        # fare taramasında kaçırmak başarısızlık değildir.
        "host_coverage": host_coverage(),
        "scope_note": SCOPE_NOTE,
        "missing_sequences": dict(MISSING_SEQUENCES),
    }


def verify_coords() -> list[str]:
    """coords alanlarını LipL32 referans dizisine karşı doğrular.

    Dönen liste BOŞ ise her şey tutarlı. Değilse her satır bir uyuşmazlıktır.
    Bu, 'koordinat uydurma' hatasına karşı makine kontrolüdür.
    """
    problems: list[str] = []
    ref = lipl32_seq()
    if not ref:
        return ["refdata/leptospira/O34094_LipL32.fasta yok — koordinat "
                "doğrulanamadı (dürüst 'çalıştırılmadı')"]
    for p in PEPTIDES:
        if p.coords is None or p.uniprot != "O34094":
            continue
        start, end = (int(x) for x in p.coords.split("-"))
        got = ref[start - 1:end]
        if got != p.seq:
            problems.append(
                f"{p.name}: coords {p.coords} -> '{got}' ama dizi '{p.seq}'"
            )
        if end - start + 1 != len(p.seq):
            problems.append(
                f"{p.name}: coords uzunluğu {end - start + 1} != dizi {len(p.seq)}"
            )
    return problems


def find_overlaps(seq: str, k: int = 8) -> list[tuple[ValPeptide, str]]:
    """Bir aday peptidi doğrulama setiyle eşleştirir.

    iedb_match.py ile aynı kanıt sınıflandırmasını kullanır:
      exact / aday_epitop_içinde / epitop_aday_içinde -> 'validated'
      yalnız ortak k-mer çekirdek                      -> 'weak'
    Dönen: [(ValPeptide, match_type), ...]
    """
    s = (seq or "").strip().upper()
    hits: list[tuple[ValPeptide, str]] = []
    if not s:
        return hits
    for p in PEPTIDES:
        e = p.seq.upper()
        if s == e:
            hits.append((p, "exact"))
        elif e in s:
            hits.append((p, "epitop_aday_içinde"))
        elif s in e:
            hits.append((p, "aday_epitop_içinde"))
        else:
            kmers = {s[i:i + k] for i in range(len(s) - k + 1)} if len(s) >= k else set()
            if kmers and any(e[i:i + k] in kmers for i in range(len(e) - k + 1)):
                hits.append((p, f"ortak_{k}mer"))
    return hits
