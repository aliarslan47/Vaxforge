"""İmmün korelat uyarıları — pipeline'ın NE ÖLÇMEDİĞİNİ raporda zorunlu kılar.

NEDEN BU MODÜL VAR
------------------
Bu pipeline epitop/protein düzeyinde ÖNGÖRÜ üretir. Bir aşının işe yarayıp
yaramadığını belirleyen şey ise o patojende KORUMANIN KORELATIDIR. İkisi her
organizmada aynı şey değildir ve bazı organizmalarda korelat, bu pipeline'ın
yapısal olarak göremeyeceği bir moleküle (ör. LPS gibi bir karbonhidrata)
bağlıdır.

Bu repodaki araç-dürüstlüğü deseninin yeni bir varyantı şudur: ölçtüğümüz
nicelik (immünojenisite / antikor tanıması) hedef niceliğin (koruma) yalnız
zayıf vekili değil, BAZI DURUMLARDA TERS İŞARETLİDİR. Örnek: ters vaksinolojiyle
bulunan bir Leptospira antijeni güçlü IgG uyardı, korumadı ve aşılanan grup
klinik belirtileri kontrolden DAHA HIZLI geliştirdi (Barbosa 2024).

Bu yüzden korelat uyarıları bir "ek not" değil, raporun ZORUNLU bölümüdür:
caveats() boş dönmüyorsa rapor o bölümü yazmak zorundadır.

TASARIM
-------
  * Uyarılar organizma profiline göre kayıtlıdır; profil tanınmıyorsa yalnız
    her organizma için geçerli GENEL uyarılar döner (sessiz kalmaz).
  * Her uyarı atıflıdır (DOI). Atıfsız uyarı yazılmaz.
  * `blocks_claim` alanı, o uyarının raporda hangi iddiayı GEÇERSİZ kıldığını
    makine-okunur biçimde söyler; rapor bu iddiayı kurmamalıdır.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

METHOD = "GERÇEK (hakemli literatür, korelat kanıtı)"


@dataclass(frozen=True)
class Caveat:
    """Tek bir korelat uyarısı."""

    id: str
    severity: str            # "critical" | "major" | "info"
    title_tr: str
    title_en: str
    text_tr: str
    text_en: str
    citations: tuple[str, ...] = field(default_factory=tuple)
    blocks_claim: str = ""   # bu uyarı hangi iddiayı geçersiz kılar

    def title(self, lang: str = "tr") -> str:
        return self.title_en if lang == "en" else self.title_tr

    def text(self, lang: str = "tr") -> str:
        return self.text_en if lang == "en" else self.text_tr


# --------------------------------------------------------------------------- #
# Her organizma için geçerli genel uyarılar
# --------------------------------------------------------------------------- #
_GENERAL: list[Caveat] = [
    Caveat(
        id="prediction_is_not_protection",
        severity="info",
        title_tr="Öngörü koruma değildir",
        title_en="Prediction is not protection",
        text_tr=(
            "Bu rapordaki sıralama, epitop/antijen ÖNGÖRÜSÜDÜR. Hiçbir satır "
            "koruma iddiası taşımaz; koruma yalnız uygun konakta challenge "
            "deneyiyle gösterilebilir. Rapor, test edilecek hipotezleri "
            "önceliklendirir."
        ),
        text_en=(
            "The ranking in this report is an epitope/antigen PREDICTION. No row "
            "claims protection; protection can only be shown by a challenge "
            "experiment in a suitable host. The report prioritises hypotheses "
            "to be tested."
        ),
        blocks_claim="Bu aday koruyucudur.",
    ),
]

# --------------------------------------------------------------------------- #
# Leptospira — taramada belgelenen üç kritik uyarı
# (ayrıntılı gerekçe: leptospira literatür taraması, 2026-11-22)
# --------------------------------------------------------------------------- #
_LEPTOSPIRA: list[Caveat] = [
    Caveat(
        id="lepto_lps_correlate_unmodelled",
        severity="critical",
        title_tr="Doğrulanmış korelat anti-LPS'tir; bu pipeline onu modellemez",
        title_en="The validated correlate is anti-LPS; this pipeline does not model it",
        text_tr=(
            "Leptospira'da deneysel olarak kanıtlanmış ve DÜZENLEYİCİ OLARAK "
            "KABUL EDİLMİŞ tek koruma korelatı, LPS'ye özgü AGLÜTİNE EDEN "
            "IgG'dir. Pasif bağışıklamada aglütine eden monoklonal antikorlar "
            "(22 kDa L-LPS) kobayları ölümcül pulmoner hemorajiden tamamen "
            "korudu; aglütine ETMEYEN monoklonal (8 kDa) korumadı. Bu korelat, "
            "hamster challenge potency testinin yerine geçen USDA ve Ph. Eur. "
            "ELISA'larının temelidir. Söz konusu korelat bir KARBONHİDRATA "
            "dayanır; protein/epitop tabanlı bu pipeline onu yapısal olarak "
            "göremez. Dolayısıyla bu rapordaki protein adayları, bilinen ana "
            "koruma mekanizmasının DIŞINDAKİ bir mekanizmayı hedefler."
        ),
        text_en=(
            "In Leptospira the only experimentally proven and REGULATORY-ACCEPTED "
            "correlate of protection is AGGLUTINATING IgG specific for LPS. In "
            "passive immunisation, agglutinating monoclonals (22 kDa L-LPS) fully "
            "protected guinea pigs from fatal pulmonary haemorrhage, whereas a "
            "non-agglutinating monoclonal (8 kDa) did not. This correlate "
            "underpins the USDA and Ph. Eur. ELISAs that replaced the hamster "
            "challenge potency test. That correlate rests on a CARBOHYDRATE; a "
            "protein/epitope-based pipeline cannot represent it. The protein "
            "candidates here therefore target a mechanism OUTSIDE the known "
            "principal route to protection."
        ),
        citations=(
            "Challa S, ve ark. Vaccine. 2011;29(27):4431-4. doi:10.1016/j.vaccine.2011.04.041",
            "Yan KT, ve ark. Vet Microbiol. 1999;69(3):173-87. doi:10.1016/s0378-1135(99)00103-0",
            "Ruby KW, Srinivas GB. Biologicals. 2013;41(5):308-14. doi:10.1016/j.biologicals.2013.06.007",
            "Klaasen HLBM, ve ark. Biologicals. 2013;41(5):315-22. doi:10.1016/j.biologicals.2013.06.012",
            "Adler B. Curr Top Microbiol Immunol. 2015;387:251-72. doi:10.1007/978-3-662-45059-8_10",
        ),
        blocks_claim="Protein adayı bilinen koruma korelatını hedefliyor.",
    ),
    Caveat(
        id="lepto_titre_does_not_predict",
        severity="critical",
        title_tr="Antikor titresi korumayı öngörmez",
        title_en="Antibody titre does not predict protection",
        text_tr=(
            "Leptospira'da antikor düzeyi ile koruma arasındaki kopukluk birden "
            "çok platformda ve bağımsız gruplarda gösterilmiştir. BCG-vektörlü "
            "aşılar hamsterları %80-100 korudu ve bazı gruplarda renal "
            "kolonizasyonu da engelledi; buna karşın yazarlar hümoral yanıtın "
            "koruma ile İLİŞKİLİ OLMADIĞINI bildirdi. İlk mRNA-LNP leptospiroz "
            "aşısı, rekombinant proteinden daha yüksek IgG, in vitro kompleman "
            "aracılı öldürme ve IFN-γ ELISpot ile Th1 yanıtı üretmesine rağmen "
            "hamsterda en iyi %25 sağkalım sağladı. Köpek ticari aşılarının "
            "sistematik derlemesi de aşıyla uyarılan MAT titrelerinin süresinin "
            "koruma ile korele olmadığını bildiriyor. SONUÇ: antikor "
            "tanıması/titresi üzerinden kurulan bir etiket, bu başarılı aşıları "
            "BAŞARISIZ olarak işaretlerdi."
        ),
        text_en=(
            "In Leptospira the disconnect between antibody level and protection "
            "has been shown across platforms and by independent groups. "
            "BCG-vectored vaccines protected 80-100% of hamsters, some also "
            "preventing renal colonisation, yet the authors reported that the "
            "humoral response was NOT related to protection. The first mRNA-LNP "
            "leptospirosis vaccine produced higher IgG than recombinant protein, "
            "in vitro complement-mediated killing and a Th1 IFN-γ ELISpot "
            "response, yet reached at best 25% survival in hamsters. A systematic "
            "review of commercial canine vaccines likewise reports that the "
            "duration of vaccine-induced MAT titres does not correlate with "
            "protection. CONSEQUENCE: a label built on antibody "
            "recognition/titre would mark these successful vaccines as FAILURES."
        ),
        citations=(
            "Oliveira TL, ve ark. Vaccine. 2019;37(6):776-82. doi:10.1016/j.vaccine.2018.12.050",
            "Dorneles J, ve ark. Vaccine. 2020;38(51):8136-44. doi:10.1016/j.vaccine.2020.10.086",
            "Techawiwattanaboon T, ve ark. Vaccine. 2026;73:128099. doi:10.1016/j.vaccine.2025.128099",
            "Bergmann Esteves S, ve ark. Vaccine. 2022;40(12):1722-40. doi:10.1016/j.vaccine.2022.02.021",
        ),
        blocks_claim="Yüksek antikor yanıtı bu adayın koruyucu olacağını gösterir.",
    ),
    Caveat(
        id="lepto_immunogenicity_can_harm",
        severity="critical",
        title_tr="Yüksek immünojenisite hastalık ağırlaşmasıyla birlikte görülebilir",
        title_en="High immunogenicity can come with disease enhancement",
        text_tr=(
            "Ters vaksinolojiyle bulunan bir Leptospira metil-kabul eden "
            "kemotaksi proteini (MCP), C3H/HeJ farelerinde tek dozla anlamlı "
            "IgG yanıtı (P<0,0001; özellikle IgG1 ve IgG2a) uyardı. Buna karşın "
            "ne hastalıktan ne renal kolonizasyondan korudu; dahası AŞILANAN "
            "GRUP klinik belirtileri kontrol grubundan DAHA HIZLI geliştirdi, "
            "sağkalım eğrileri anlamlı biçimde farklıydı ve makroskopik/"
            "histopatolojik lezyonlar değişti. SONUÇ: 'yüksek immünojenisite = "
            "iyi aday' kuralı bu organizmada güvenli değildir; bu rapordaki "
            "immünojenisite skorları bir GÜVENLİK güvencesi değildir."
        ),
        text_en=(
            "A Leptospira methyl-accepting chemotaxis protein (MCP) identified by "
            "reverse vaccinology elicited a significant IgG response in C3H/HeJ "
            "mice after a single dose (P<0.0001; notably IgG1 and IgG2a). Yet it "
            "protected against neither disease nor renal colonisation, and the "
            "VACCINATED group developed clinical signs FASTER than controls, with "
            "significantly different survival curves and altered gross and "
            "histopathological lesions. CONSEQUENCE: the rule 'high "
            "immunogenicity = good candidate' is not safe in this organism; the "
            "immunogenicity scores in this report are not a SAFETY assurance."
        ),
        citations=(
            "Barbosa LN, ve ark. PLoS Negl Trop Dis. 2024;18(9):e0012155. "
            "doi:10.1371/journal.pntd.0012155",
        ),
        blocks_claim="İmmünojenisitesi yüksek aday güvenli/iyi adaydır.",
    ),
    Caveat(
        id="lepto_surface_label_unreliable",
        severity="major",
        title_tr="Yüzey maruziyeti etiketi bu türde yöntem-bağımlıdır",
        title_en="The surface-exposure label is method-dependent in this organism",
        text_tr=(
            "Leptospira membran fraksiyonlama yöntemleri belirli dış membran "
            "protein sınıflarına (lipoprotein ile transmembran β-varil) seçicidir "
            "ve YANLIŞ SINIFLANDIRMAYA yol açabilir; spiroket dış membranı "
            "kırılgandır ve dış membranda bulunmak yüzeyde maruz olmayı garanti "
            "etmez. Ayrıca antikorlar genellikle rekombinant proteine karşı "
            "üretildiği için native yüzey epitoplarını tanımayabilir. Bu nedenle "
            "yüzey/yerleşim öngörüsü bu organizmada SERT bir eleme ölçütü olarak "
            "kullanılmamalıdır. Bu repodaki benchmark'ta da LigA, PSORTb'nin "
            "'cytoplasm' ataması (skor 9,97) nedeniyle kaçırılmıştır — oysa "
            "literatürde yerleşik bir yüzey lipoproteinidir."
        ),
        text_en=(
            "Methods for fractionating Leptospira membranes are selective for "
            "particular classes of outer-membrane protein (lipoproteins versus "
            "transmembrane β-barrels) and can lead to MISCLASSIFICATION; the "
            "spirochaetal outer membrane is fragile, and presence in the outer "
            "membrane does not guarantee surface exposure. Antibodies are usually "
            "raised against recombinant protein and may not recognise native "
            "surface epitopes. Surface/localisation prediction should therefore "
            "not be used as a HARD filter in this organism. In this repository's "
            "own benchmark, LigA was missed because PSORTb assigned it "
            "'cytoplasm' (score 9.97), although it is an established surface "
            "lipoprotein in the literature."
        ),
        citations=(
            "Pinne M, Haake D. J Vis Exp. 2011;(53):2805. doi:10.3791/2805",
        ),
        blocks_claim="Yüzeyde öngörülmeyen protein aday değildir.",
    ),
    Caveat(
        id="lepto_abundance_conservation_misleads",
        severity="major",
        title_tr="Bolluk + korunmuşluk + immünodominans koruma öngörmez (LipL32)",
        title_en="Abundance + conservation + immunodominance do not predict protection (LipL32)",
        text_tr=(
            "LipL32 dış membran proteininin %75'ine kadarını oluşturur, yalnız "
            "patojenik türlerde bulunur, yüksek korunmuştur ve doğal enfeksiyonda "
            "baskın immünojendir — yani bir ters vaksinoloji sıralamasının en üst "
            "sırasına oturacak bütün özelliklere sahiptir. Buna karşın transpozon "
            "lipL32 mutantı hamsterda akut hastalıkta ve sıçanda kronik "
            "kolonizasyonda yabanıl tipten AYIRT EDİLEMEZ; geniş adjuvan ve "
            "vektör yelpazesine rağmen aşı sonuçları karışıktır. Bu rapordaki "
            "bolluk/korunmuşluk/antijenisite sıralaması bu organizmada yanlış "
            "adayı ilk sıraya koyabilir."
        ),
        text_en=(
            "LipL32 accounts for up to 75% of outer-membrane protein, occurs only "
            "in pathogenic species, is highly conserved and is the dominant "
            "immunogen in natural infection — every property that would place it "
            "at the top of a reverse-vaccinology ranking. Yet a transposon lipL32 "
            "mutant is INDISTINGUISHABLE from wild type in acute hamster disease "
            "and chronic rat colonisation, and vaccine results remain mixed across "
            "a wide range of adjuvants and vectors. The abundance/conservation/"
            "antigenicity ranking in this report may therefore place the wrong "
            "candidate first in this organism."
        ),
        citations=(
            "Murray GL, ve ark. Infect Immun. 2009;77(3):952-8. doi:10.1128/IAI.01370-08",
            "Murray GL. Vet Microbiol. 2013;162(2-4):305-14. doi:10.1016/j.vetmic.2012.11.005",
        ),
        blocks_claim="En bol ve en korunmuş yüzey proteini en iyi adaydır.",
    ),
    Caveat(
        id="lepto_published_protection_noisy",
        severity="major",
        title_tr="Yayınlanmış koruma etiketleri metodolojik olarak gürültülü",
        title_en="Published protection labels are methodologically noisy",
        text_tr=(
            "Alanın kendi değerlendirmesine göre, YETERSİZ CHALLENGE DOZLARI ve "
            "UYGUNSUZ İSTATİSTİKSEL ANALİZ kullanılması, rekombinant proteinlerle "
            "elde edilen koruma iddialarının çoğunu geçersiz kılar. Bu nedenle "
            "literatürden toplanan 'korudu/korumadı' etiketleri kalite ağırlığı "
            "(challenge dozu LD50 cinsinden, grup büyüklüğü, istatistiksel test, "
            "homolog/heterolog challenge) olmadan kullanılmamalıdır."
        ),
        text_en=(
            "By the field's own assessment, the use of INADEQUATE CHALLENGE DOSES "
            "and INAPPROPRIATE STATISTICAL ANALYSIS invalidates many of the "
            "claims of protection obtained with recombinant proteins. Labels of "
            "'protected / not protected' harvested from the literature should "
            "therefore not be used without a quality weight (challenge dose in "
            "LD50, group size, statistical test, homologous versus heterologous "
            "challenge)."
        ),
        citations=(
            "Adler B. Curr Top Microbiol Immunol. 2015;387:251-72. "
            "doi:10.1007/978-3-662-45059-8_10",
        ),
        blocks_claim="Literatürdeki koruma etiketleri doğrudan eğitim verisi olarak kullanılabilir.",
    ),
    Caveat(
        id="lepto_iedb_empty",
        severity="major",
        title_tr="IEDB bu organizma için pratikte boştur — recall ölçülemez",
        title_en="IEDB is practically empty for this organism — recall cannot be measured",
        text_tr=(
            "IEDB IQ-API 2026-11-22'de sorgulandı: Leptospira cins düğümünde 2, "
            "L. interrogans'ta 14 (11'i non-negatif lineer), L. kirschneri'de 1, "
            "diğer patojenik türlerde 0 yapı kaydı. Toplam ~16. Bu kayıtların "
            "çoğu tek bir TANISAL fay gösterimi çalışmasından gelir ve GroEL "
            "kaydının kaynağı bir insan HSP60 otoantikor çalışmasıdır. Bu nedenle "
            "IEDB tabanlı recall benchmark'ı bu organizmada ANLAMSIZDIR; "
            "lepto_validation.py modülündeki literatür temelli set kullanılır ve "
            "payda raporda açıkça yazılır."
        ),
        text_en=(
            "The IEDB IQ-API was queried on 2026-11-22: 2 structure records at "
            "the Leptospira genus node, 14 for L. interrogans (11 non-negative "
            "linear), 1 for L. kirschneri and 0 for the other pathogenic species "
            "— about 16 in total. Most come from a single DIAGNOSTIC phage-display "
            "study, and the GroEL record originates from a human HSP60 "
            "autoantibody study. An IEDB-based recall benchmark is therefore "
            "MEANINGLESS here; the literature-based set in lepto_validation.py is "
            "used instead and its denominator is stated explicitly in the report."
        ),
        citations=(
            "Vita R, ve ark. The Immune Epitope Database (IEDB): 2024 update. "
            "Nucleic Acids Res. 2025;53(D1):D436-D443. doi:10.1093/nar/gkae1092",
            "Sayımlar literatürden alınmadı, IEDB IQ-API'den 2026-11-22'de CANLI "
            "ÖLÇÜLDÜ: GET https://query-api.iedb.org/epitope_search?"
            "structure_type=eq.Linear+peptide&qualitative_measure=neq.Negative&"
            "source_organism_iris=cs.{NCBITaxon:<id>} (Prefer: count=exact, "
            "Content-Range başlığı).",
        ),
        blocks_claim="IEDB recall oranı bu organizmada pipeline başarımını ölçer.",
    ),
    Caveat(
        id="lepto_ifn_host_gap",
        severity="info",
        title_tr="IFN-γ öngörüsü ilgili konaklar için modelsizdir",
        title_en="IFN-γ prediction has no model for the relevant hosts",
        text_tr=(
            "Leptospira'da alanın en iyi doğrulanmış PROTEİN korelatı sığır "
            "Hardjo modelindeki Th1/IFN-γ yanıtıdır (IFN-γ üreten CD4+ ve WC1+ "
            "γδ T hücreleri; koruyan monovalan aşı bu yanıtı uyarır, korumayan "
            "pentavalan uyarmaz). Ancak kullanılan IFN-γ indükleme yordayıcısının "
            "modeli yalnız İNSAN ve FARE içindir; sığır ve hamster için model "
            "YOKTUR. Bu nedenle bu konaklarda IFN-γ alanı 'ölçülmedi' olarak "
            "raporlanır ve bu bir NEGATİF BULGU DEĞİLDİR."
        ),
        text_en=(
            "The best-validated PROTEIN correlate in the Leptospira field is the "
            "Th1/IFN-γ response in the bovine Hardjo model (IFN-γ-producing CD4+ "
            "and WC1+ γδ T cells; the protective monovalent vaccine induces it, "
            "the non-protective pentavalent one does not). However, the IFN-γ "
            "induction predictor used here has models only for HUMAN and MOUSE; "
            "there is NO model for cattle or hamster. The IFN-γ field is therefore "
            "reported as 'not measured' for those hosts, and that is NOT a "
            "negative finding."
        ),
        citations=(
            "Naiman BM, ve ark. Infect Immun. 2001;69(12):7550-8. "
            "doi:10.1128/IAI.69.12.7550-7558.2001",
            "Brown RA, ve ark. Vaccine. 2003;21(27-30):4448-58. "
            "doi:10.1016/s0264-410x(03)00439-0",
            "Deveson Lucas DS, ve ark. Vet Microbiol. 2014;169(3-4):163-70. "
            "doi:10.1016/j.vetmic.2013.12.016",
        ),
        blocks_claim="IFN-γ ölçülmemiş olması adayın Th1 uyarmadığını gösterir.",
    ),
]

REGISTRY: dict[str, list[Caveat]] = {
    "leptospira": _LEPTOSPIRA,
}

# Profil adı eşleştirme: kullanıcı 'leptospira', 'Leptospira interrogans',
# 'lepto' vb. yazabilir. Anahtar kelime -> kayıt anahtarı.
_ALIASES = {
    "lepto": "leptospira",
    "leptospira": "leptospira",
    "leptospirosis": "leptospira",
    "leptospiroz": "leptospira",
    "interrogans": "leptospira",
    "borgpetersenii": "leptospira",
    "kirschneri": "leptospira",
}


def resolve(profile: str | None) -> str | None:
    """Serbest metin profil/organizma adını kayıt anahtarına çözer."""
    if not profile:
        return None
    s = str(profile).strip().lower()
    for token, key in _ALIASES.items():
        if token in s:
            return key
    return None


def detect_from_input(path: str | None = None, extra: str | None = None,
                      max_bytes: int = 200_000) -> str | None:
    """Girdi dosyasının adı + başlıklarından organizmayı çıkarır.

    Profil seçimi ('bacteria', 'virus') organizmayı söylemez; korelat uyarıları
    ORGANİZMAYA bağlıdır. Bu yüzden organizma, dosya adı ve FASTA/GenBank
    başlıklarından (ORGANISM/DEFINITION/>...) okunur. Bulunamazsa None döner ve
    yalnız genel uyarılar gösterilir — sessiz kalınmaz.
    """
    blobs: list[str] = []
    if extra:
        blobs.append(str(extra))
    if path:
        blobs.append(str(path))
        try:
            p = Path(path)
            if p.is_file():
                if p.suffix == ".gz":
                    import gzip
                    with gzip.open(p, "rt", errors="replace") as fh:
                        blobs.append(fh.read(max_bytes))
                else:
                    with open(p, errors="replace") as fh:
                        blobs.append(fh.read(max_bytes))
        except OSError:
            pass
    for b in blobs:
        key = resolve(b)
        if key:
            return key
    return None


def caveats(profile: str | None, include_general: bool = True) -> list[Caveat]:
    """Bu organizma için geçerli korelat uyarıları, önem sırasına göre.

    Profil tanınmasa bile GENEL uyarılar döner — modül hiçbir koşulda sessiz
    kalmaz.
    """
    key = resolve(profile)
    out = list(REGISTRY.get(key, [])) if key else []
    if include_general:
        out = out + list(_GENERAL)
    order = {"critical": 0, "major": 1, "info": 2}
    return sorted(out, key=lambda c: order.get(c.severity, 9))


def is_mandatory(profile: str | None) -> bool:
    """Bu organizmada korelat bölümü ZORUNLU mu (critical uyarı var mı)?"""
    return any(c.severity == "critical" for c in caveats(profile))


def summary(profile: str | None) -> dict:
    cs = caveats(profile)
    return {
        "profile": profile,
        "resolved": resolve(profile),
        "method": METHOD,
        "mandatory": is_mandatory(profile),
        "n_caveats": len(cs),
        "n_critical": len([c for c in cs if c.severity == "critical"]),
        "n_major": len([c for c in cs if c.severity == "major"]),
        "ids": [c.id for c in cs],
        "blocked_claims": [c.blocks_claim for c in cs if c.blocks_claim],
    }
