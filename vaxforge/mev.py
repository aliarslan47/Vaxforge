"""Multi-epitop aşı (MEV) konstruktu inşası ve karakterizasyonu.

Sıralanmış tekil epitopları (B / MHC-I / MHC-II) standart linker'lar ve bir
adjuvan ile tek bir protein dizisine birleştirir; ardından konstrukt-seviyesi
fizikokimyasal + immünoinformatik özellikleri hesaplar.

Linker kuralları (literatür standardı):
  - Adjuvan → ilk epitop : EAAAK  (rijit, adjuvanı ayrı katlar)
  - MHC-I (CTL) arası      : AAY    (proteazom kesim bölgesi)
  - MHC-II (HTL) arası     : GPGPG  (HTL işlemeyi kolaylaştırır)
  - B-hücre (LBL) arası    : KK     (bağımsız immünojenite)
  - Bloklar arası          : GPGPG

Karakterizasyon mevcut predictor'ları yeniden kullanır (protparam, antigen_acc,
allergen, toxinpred, funnel homoloji). Solubility / 2° yapı / disorder BİLEREK
kapsam dışıdır (ayrı predictor gerektirir).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import (allergen, antigen_acc, disorder, iapred, secstruct, sequtils,
               solubility, toxinpred)
from .models import Peptide, ProteinRecord

# ── Adjuvan kütüphanesi (N-ucuna EAAAK ile eklenir) ────────────────────────
# Her dizi literatürden/UniProt'tan doğrulanmıştır (uydurma YOK — VaxForge kuralı).
# Kullanıcı hangisini kullanacağını seçebilir; ADJUVANT_INFO metaverisi arayüze
# ve rapora aktarılır (atıf her zaman yanında).
ADJUVANTS: dict[str, str] = {
    # İnsan β-defensin (hBD) — TLR4 agonisti, DC olgunlaşması; en yaygın MEV adjuvanı
    "beta_defensin": "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK",
    # 50S ribozomal protein L7/L12 (M. tuberculosis) — TLR4 agonisti
    "l7_l12": (
        "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVEAAE"
        "EQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEA"
        "KAKLEAAGATVTVK"
    ),
    # PADRE — pan-HLA-DR bağlayan universal CD4+ Th epitopu (Alexander 1994)
    "padre": "AKFVAAWTLKAAA",
    # RS09 — sentetik TLR4 agonist peptid (Shanmugam 2012)
    "rs09": "APPHALS",
    # HBHA — heparin-bağlayan hemaglutinin adezin (M. tuberculosis, UniProt P9WIP9);
    # TLR4 üzerinden T-hücre aktivasyonu ve IFN-γ (Menozzi 1996; Pethe 2001)
    "hbha": (
        "MAENSNIDDIKAPLLAALGAADLALATVNELITNLRERAEETRTDTRSRVEESRARLTK"
        "LLQEDLPEQLTELREKFTAEELRKAAEGYLEAATSRYNELVERGEAALERLRSQQSFEE"
        "EVSARAEGYVDQAVELTQEALGTVASQTRAVGERAAKLVGIELPKKAAPAKKAAPAKKA"
        "APAKKAAPAKKAAAKKAPAKKAAAKKVTQK"
    ),
    "none": "",
}

# Adjuvan metaverisi — arayüz seçici + rapor için (label TR/EN, mekanizma, atıf anahtarı).
ADJUVANT_INFO: dict[str, dict] = {
    "beta_defensin": {"label_tr": "β-defensin (hBD)", "label_en": "β-defensin (hBD)",
                      "tlr": "TLR4", "citation": "Biragyn et al. 2002",
                      "desc_tr": "İnsan β-defensin; TLR4 üzerinden dendritik hücre olgunlaşması. En yaygın MEV adjuvanı.",
                      "desc_en": "Human β-defensin; TLR4-driven dendritic-cell maturation. Most common MEV adjuvant."},
    "l7_l12": {"label_tr": "50S ribozomal L7/L12", "label_en": "50S ribosomal L7/L12",
               "tlr": "TLR4", "citation": "Rahmani et al. 2019",
               "desc_tr": "M. tuberculosis 50S ribozomal protein L7/L12; TLR4 agonisti.",
               "desc_en": "M. tuberculosis 50S ribosomal protein L7/L12; TLR4 agonist."},
    "padre": {"label_tr": "PADRE (universal Th)", "label_en": "PADRE (universal Th)",
              "tlr": "—", "citation": "Alexander et al. 1994",
              "desc_tr": "Pan-HLA-DR bağlayan universal CD4+ yardımcı-T epitopu; CTL yanıtını güçlendirir.",
              "desc_en": "Pan-HLA-DR-binding universal CD4+ helper-T epitope; boosts CTL response."},
    "rs09": {"label_tr": "RS09 (sentetik TLR4)", "label_en": "RS09 (synthetic TLR4)",
             "tlr": "TLR4", "citation": "Shanmugam et al. 2012",
             "desc_tr": "Sentetik TLR4 agonist heptapeptid (APPHALS); LPS yerine güvenli mimik.",
             "desc_en": "Synthetic TLR4-agonist heptapeptide (APPHALS); safe LPS mimic."},
    "hbha": {"label_tr": "HBHA (M. tuberculosis)", "label_en": "HBHA (M. tuberculosis)",
             "tlr": "TLR4", "citation": "Pethe et al. 2001",
             "desc_tr": "Heparin-bağlayan hemaglutinin adezin; T-hücre aktivasyonu ve IFN-γ salımı.",
             "desc_en": "Heparin-binding hemagglutinin adhesin; T-cell activation and IFN-γ release."},
    "none": {"label_tr": "Adjuvan yok", "label_en": "No adjuvant",
             "tlr": "—", "citation": "",
             "desc_tr": "Adjuvansız — yalnız epitop kaseti (karşılaştırma / dış adjuvanla formülasyon için).",
             "desc_en": "No adjuvant — epitope cassette only (for comparison / external-adjuvant formulation)."},
}


def list_adjuvants() -> list[dict]:
    """Arayüz seçici için adjuvan kataloğu (anahtar + metaveri + uzunluk)."""
    out = []
    for key, seq in ADJUVANTS.items():
        info = ADJUVANT_INFO.get(key, {})
        out.append({
            "key": key,
            "label_tr": info.get("label_tr", key),
            "label_en": info.get("label_en", key),
            "tlr": info.get("tlr", "—"),
            "citation": info.get("citation", ""),
            "desc_tr": info.get("desc_tr", ""),
            "desc_en": info.get("desc_en", ""),
            "length": len(seq),
        })
    return out

# ── Linker'lar ─────────────────────────────────────────────────────────────
LINKERS = {
    "adjuvant": "EAAAK",
    "MHC-I": "AAY",
    "MHC-II": "GPGPG",
    "B": "KK",
    "block": "GPGPG",
}


@dataclass
class MEVConstruct:
    """İnşa edilen konstrukt + izlenebilir kompozisyon."""

    seq: str
    adjuvant: str                      # kullanılan adjuvan anahtarı
    cassette: str = ""                 # adjuvan HARİÇ epitop+linker kaseti
    components: list[dict] = field(default_factory=list)   # sıralı parça izi
    n_by_kind: dict[str, int] = field(default_factory=dict)


def _pick(peptides: list[Peptide], kind: str, top: int) -> list[Peptide]:
    """Verilen tipte, filtreden geçmiş, candidacy'ye göre en iyi `top` epitop."""
    pool = [p for p in peptides if p.kind == kind and p.passed]
    pool.sort(key=lambda p: p.candidacy, reverse=True)
    # aynı diziyi iki kez koymayalım
    seen: set[str] = set()
    out: list[Peptide] = []
    for p in pool:
        s = sequtils.sanitize(p.seq)
        if s and s not in seen:
            seen.add(s)
            out.append(p)
        if len(out) >= top:
            break
    return out


def build_construct(
    peptides: list[Peptide],
    adjuvant: str = "beta_defensin",
    top_i: int = 6,
    top_ii: int = 4,
    top_b: int = 4,
) -> MEVConstruct:
    """Top epitopları linker + adjuvan ile birleştirip konstrukt üretir.

    Sıra:  [adjuvan]-EAAAK- CTL(AAY) -GPGPG- HTL(GPGPG) -GPGPG- LBL(KK)
    """
    adj_seq = ADJUVANTS.get(adjuvant, ADJUVANTS["beta_defensin"])
    ctl = _pick(peptides, "MHC-I", top_i)
    htl = _pick(peptides, "MHC-II", top_ii)
    lbl = _pick(peptides, "B", top_b)

    parts: list[str] = []
    comp: list[dict] = []

    def add(seq: str, role: str, src: str = "") -> None:
        parts.append(seq)
        comp.append({"seq": seq, "role": role, "source": src})

    if adj_seq:
        add(adj_seq, "adjuvant", adjuvant)
        add(LINKERS["adjuvant"], "linker")

    for block, peps, lk in (
        ("MHC-I", ctl, LINKERS["MHC-I"]),
        ("MHC-II", htl, LINKERS["MHC-II"]),
        ("B", lbl, LINKERS["B"]),
    ):
        if not peps:
            continue
        if parts and comp[-1]["role"] != "linker":
            add(LINKERS["block"], "linker")
        for i, p in enumerate(peps):
            if i > 0:
                add(lk, "linker")
            add(sequtils.sanitize(p.seq), f"epitope:{block}", p.parent)

    seq = "".join(parts)
    # Kaset = adjuvan ve onu izleyen tek EAAAK ayırıcı HARİÇ kalan her şey.
    # İnsan-homoloji kontrolü bu kaset üzerinde yapılır: adjuvan (β-defensin)
    # bilerek insan-kökenlidir, dolayısıyla self-benzerliği yorumu yanıltmasın.
    cassette_parts, dropped_eaaak = [], False
    for c in comp:
        if c["role"] == "adjuvant":
            continue
        if not dropped_eaaak and c["role"] == "linker" and c["seq"] == LINKERS["adjuvant"]:
            dropped_eaaak = True   # adjuvanı ayıran tek EAAAK'ı da at
            continue
        cassette_parts.append(c["seq"])
    cassette = "".join(cassette_parts)
    return MEVConstruct(
        seq=seq,
        adjuvant=adjuvant,
        cassette=cassette,
        components=comp,
        n_by_kind={"MHC-I": len(ctl), "MHC-II": len(htl), "B": len(lbl)},
    )


# ── Fizikokimyasal yardımcılar ─────────────────────────────────────────────
def aliphatic_index(seq: str) -> float:
    """Alifatik indeks (Ikai 1980) — termostabilite göstergesi.

    AI = X_Ala + 2.9·X_Val + 3.9·(X_Ile + X_Leu)   (X = mol %, 0-100)
    """
    s = sequtils.sanitize(seq)
    if not s:
        return 0.0
    n = len(s)
    a = 100 * s.count("A") / n
    v = 100 * s.count("V") / n
    i = 100 * s.count("I") / n
    ll = 100 * s.count("L") / n
    return round(a + 2.9 * v + 3.9 * (i + ll), 2)


def _charged(seq: str) -> tuple[int, int]:
    """(negatif Asp+Glu, pozitif Arg+Lys) rezidü sayıları."""
    s = sequtils.sanitize(seq)
    neg = s.count("D") + s.count("E")
    pos = s.count("R") + s.count("K")
    return neg, pos


def _human_similarity(seq: str) -> dict:
    """Konstruktun insan proteomuna en iyi %kimliği (diamond, funnel altyapısı).

    Araç/DB yoksa {'available': False}. Aksi halde en iyi %kimlik + karar.
    """
    from . import funnel  # geç import: diamond yolları modül yüklenince çözülür

    if not funnel.human_db_available():
        return {"available": False}
    rec = ProteinRecord(id="MEV", seq=seq)
    best = funnel._homology_vs_db([rec], funnel._HUMAN_DB, evalue=1e-3)
    pid = best.get("MEV", 0.0)
    return {"available": True, "best_pident": round(pid, 1),
            "similar": pid >= 35.0}


def characterize(seq: str, cassette: str | None = None) -> dict:
    """Konstrukt için özellik tablosu (kullanıcı tablosuyla aynı satırlar).

    `cassette` verilirse insan-homoloji kontrolü adjuvan hariç kaset üzerinde
    yapılır (adjuvan insan-kökenli olabilir; self-benzerlik yorumu yanıltmasın).
    """
    s = sequtils.sanitize(seq)
    pp = sequtils.protparam(s)
    neg, pos = _charged(s)

    # Antijenite — pipeline birincil aracı IApred (Miles et al. 2025);
    # yoksa ACC modeli (VaxiJen-tarzı proxy); o da yoksa fizikokimyasal proxy.
    if iapred.available():
        r = iapred.predict([("MEV", s)]).get("MEV", {})
        ag = round(r.get("norm", 0.0), 3)          # 0-1'e normalize
        ag_raw = round(r.get("score", 0.0), 3)     # ham skor (≈ -3..+3)
        ag_flag = bool(r.get("antigenic", False))  # ham > 0.3
        ag_method = f"IApred (ham={ag_raw}, {r.get('category','')})"
    elif antigen_acc.available():
        ag, ag_method = antigen_acc.predict(s), "ACC modeli (VaxiJen-tarzı proxy)"
        ag_flag = ag >= 0.5
    else:
        ag, ag_method = sequtils.antigenicity_proxy(s), "fizikokimyasal proxy"
        ag_flag = ag >= 0.4

    # Alerjenite (FAO/WHO 6-mer)
    if allergen.available():
        al = allergen.predict([s]).get(s, {})
        allergenic = bool(al.get("allergen"))
        al_method = "FAO/WHO 6-mer"
    else:
        allergenic, al_method = None, "yok"

    # Toksisite (ToxinPred2)
    if toxinpred.available():
        tx = toxinpred.predict([s]).get(s, {})
        toxic = bool(tx.get("toxic"))
        tx_method = "ToxinPred2"
    else:
        toxic, tx_method = None, "yok"

    return {
        "length": len(s),
        "mw_kda": round(pp.get("mw", 0.0) / 1000.0, 2),
        "pI": pp.get("pI"),
        "neg_residues": neg,
        "pos_residues": pos,
        "instability": pp.get("instability"),
        "stable": (pp.get("instability", 100) < 40.0),
        "aliphatic_index": aliphatic_index(s),
        "gravy": pp.get("gravy"),
        "aromaticity": pp.get("aromaticity"),
        "antigenicity": {"score": ag, "antigenic": ag_flag, "method": ag_method},
        "allergen": {"allergenic": allergenic, "method": al_method},
        "toxicity": {"toxic": toxic, "method": tx_method},
        # İçsel düzensizlik (metapredict); araç yoksa None.
        "disorder": (
            {**disorder.predict(s), "method": "metapredict V3"}
            if disorder.available() else {"percent_disordered": None, "method": "yok"}
        ),
        # İkincil yapı (S4PRED); araç yoksa None.
        "secondary_structure": (
            {**secstruct.predict(s), "method": "S4PRED"}
            if secstruct.available() else {"helix_pct": None, "method": "yok"}
        ),
        # Çözünürlük (Protein-Sol); araç yoksa None.
        "solubility": (
            {**solubility.predict(s), "method": "Protein-Sol"}
            if solubility.available() else {"scaled_solubility": None, "method": "yok"}
        ),
        # Homoloji: kaset (adjuvan hariç) üzerinde; yoksa tüm konstrukt.
        "human_similarity": _human_similarity(cassette if cassette else s),
    }


# Bu modülün kullandığı yöntemlerin citations.py'daki 'tool' etiketleri.
# Rapor bu anahtarlarla ilgili atıfları çeker (literatür daima yanında).
CITATION_KEYS = [
    "ExPASy ProtParam",              # MW, pI, GRAVY (Gasteiger 2005)
    "Kyte-Doolittle",                # GRAVY tanımı (1982)
    "Guruprasad et al. 1990",        # instability index; <40 kararlı
    "Ikai 1980",                     # aliphatic index
    "EAAAK linker (Arai et al. 2001)",       # rijit adjuvan ayırıcı
    "Livingston et al. 2002",        # GPGPG (HTL)
    "Schubert & Kohlbacher 2016",    # AAY / KK (string-of-beads)
    "metapredict V3",                # içsel düzensizlik (%disordered)
    "S4PRED (Moffat & Jones 2021)",  # ikincil yapı (%helix/strand/coil)
    "Protein-Sol (Hebditch et al. 2017)",  # çözünürlük (scaled solubility)
]


def run(
    peptides: list[Peptide],
    adjuvant: str = "beta_defensin",
    top_i: int = 6,
    top_ii: int = 4,
    top_b: int = 4,
) -> dict:
    """İnşa + karakterizasyon; rapora hazır tek sözlük döner."""
    mev = build_construct(peptides, adjuvant, top_i, top_ii, top_b)
    props = characterize(mev.seq, cassette=mev.cassette)
    # Atıflar: sabit linker/karakterizasyon + SEÇİLEN adjuvanın kendi atfı
    # (literatür daima kullanılan bileşenle eşleşir — β-defensin sabit değil).
    cites = list(CITATION_KEYS)
    adj_cite = ADJUVANT_INFO.get(mev.adjuvant, {}).get("citation")
    if adj_cite:
        cites.append(adj_cite)
    return {
        "sequence": mev.seq,
        "cassette": mev.cassette,
        "adjuvant": mev.adjuvant,
        "adjuvant_label": ADJUVANT_INFO.get(mev.adjuvant, {}).get("label_en", mev.adjuvant),
        "n_by_kind": mev.n_by_kind,
        "components": mev.components,
        "properties": props,
        "citation_keys": cites,
    }
