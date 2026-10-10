"""Adım 5 — Sağ kalım elemesi (filtre kısmı).

Sert filtre: alerjenite (non-alerjen olmalı).
Yumuşak: toksisite işaretlemesi (skor cezası), konservasyon, popülasyon kapsamı.

TOKSİSİTE NEDEN ARTIK SERT DEĞİL (2026-10-10 kalibrasyonu):
`benchmark/calibrate_filters.py` ToxinPred2'yi (Model 2 hibrit, eşik 0.6) bağımsız
kontrol setlerinde ölçtü — insan Swiss-Prot 9-mer'lerinin %75'i ve deneysel
DOĞRULANMIŞ IEDB epitoplarının %63'ü "toksik" işaretlendi; influenza M1(58-66)
GILGFVFTL bile elendi. Aynı araç TAM PROTEİNDE kusursuz çalışıyor (toksin
proteinlerinde %94 duyarlılık, non-toksinlerde %5 yanlış-pozitif). Sorun araç
değil, uzunluk sınıfı: ToxinPred2 "toxicity of *proteins*" için eğitilmiştir ve
AAC özellikleri 9-mer'de bilgi taşımaz. Sonuç: epitop aşaması rastgele bir
4'te-3 kırımına dönüşmüş, ESAT-6 / L7-L12 / ActA / LLO gibi antijenler bu yüzden
"aday peptit üretmedi" durumuna düşmüştü.

Mimari: epitop düzeyinde ToxinPred3 (peptit için eğitilmiş) yalnız İŞARETLER ve
adaylık puanını `toxicity_penalty` ile çarpar; eleme yalnız kalibre edilmiş
yüksek-güven eşiğinin (`hard_discard_score`) üstünde yapılır. Sert toksisite
elemesi araçların kendi uzunluk sınıfında, nihai MEV konstruktu üzerinde
ToxinPred2 hibrit ile korunur (bkz. mev.py).

DÜRÜSTLÜK: Gerçek araç yoksa sert filtre UYGULANMAZ. Eskiden buradaki
demonstrasyon sezgiselleri (kompozisyon tabanlı vekiller) araç eksikse sessizce
sert filtre gibi davranıp özellikle hidrofobik MHC-I bağlayıcılarını eliyordu.
Artık vekil skor yalnız ANOTASYON olarak yazılır, kapı atlanır ve 'method'
alanında açıkça belirtilir.
"""

from __future__ import annotations

from . import allergen, toxinpred3  # ToxinPred2 BİLEREK yok: protein boyunda kalibre, epitopta değil
from .config_loader import ResolvedTool
from .models import Peptide
from .sequtils import KD, sanitize


def _allergen_prob(pep: str) -> float:
    """Kompozisyon tabanlı kaba alerjen vekili — YALNIZ ANOTASYON (eleme yok)."""
    s = sanitize(pep)
    if not s:
        return 1.0
    aromatic = sum(1 for c in s if c in "FWY") / len(s)
    hydro = sum(1 for c in s if KD.get(c, 0) > 1.5) / len(s)
    return round(min(1.0, 0.3 + 0.6 * aromatic + 0.3 * hydro), 3)


def _tox_score(pep: str) -> float:
    """Sistein/pozitif yük tabanlı kaba toksisite vekili — YALNIZ ANOTASYON."""
    s = sanitize(pep)
    if not s:
        return 0.0
    cys = sum(1 for c in s if c == "C") / len(s)
    pos = sum(1 for c in s if c in "KR") / len(s)
    return round(2 * (0.6 * cys + 0.5 * pos) - 0.4, 3)


def _p(tool: ResolvedTool, name: str, default):
    rp = tool.params.get(name)
    return default if rp is None else rp.value


def run(peptides: list[Peptide], tools: dict[str, ResolvedTool]) -> tuple[list[Peptide], dict]:
    alg_tool, tox_tool = tools["allergenicity"], tools["toxicity"]
    require_non_allergen = bool(_p(alg_tool, "require_non_allergen", True))
    soft_thr = float(_p(tox_tool, "max_toxicity_score", 0.38))
    hard_thr = _p(tox_tool, "hard_discard_score", None)
    hard_thr = float(hard_thr) if hard_thr is not None else None
    penalty = float(_p(tox_tool, "toxicity_penalty", 0.9))

    active = [p for p in peptides if p.passed]
    seqs = [p.seq for p in active]

    alg_avail = allergen.available()
    alg_res = allergen.predict(seqs) if alg_avail else {}

    # Peptit düzeyi toksisite: ToxinPred3 (peptit için eğitilmiş). ToxinPred2'ye
    # DÜŞMÜYORUZ — bu uzunlukta kalibrasyonu geçmedi; yoksa kapı atlanır.
    tox_avail = toxinpred3.available()
    tox_res = toxinpred3.predict(seqs, threshold=soft_thr, model=1) if tox_avail else {}
    tox_engine = "GERÇEK (ToxinPred3 ExtraTrees, peptit düzeyi)"

    n_start = len(active)

    # --- Aşama 1: ALERJENİTE (sert — kalibrasyonda 9-mer YP ~%2-3) ----------
    after_allergen, n_allergen = [], 0
    for p in active:
        if alg_avail and p.seq in alg_res:
            is_allergen = bool(alg_res[p.seq]["allergen"])
            p.metrics["allergen"] = is_allergen
            p.metrics["allergen_match"] = alg_res[p.seq]["match"]
            p.methods["allergen"] = "GERÇEK (FAO/WHO 6-mer, UniProt allergen DB)"
        else:
            # araç yok → ELEME YOK; vekil skor yalnız anotasyon
            is_allergen = False
            p.metrics["allergen"] = None
            p.metrics["allergen_prob_proxy"] = _allergen_prob(p.seq)
            p.methods["allergen"] = "ATLANDI (gerçek alerjen DB yok; vekil yalnız anotasyon)"
        if require_non_allergen and is_allergen:
            p.passed = False
            p.notes.append("alerjen (elendi)")
            n_allergen += 1
        else:
            after_allergen.append(p)

    # --- Aşama 2: TOKSİSİTE (yumuşak işaretleme + yüksek-güven elemesi) -----
    survivors, n_toxic_flag, n_toxic_discard = [], 0, 0
    for p in after_allergen:
        if tox_avail and p.seq in tox_res:
            sc = float(tox_res[p.seq]["ml_score"])
            p.metrics["toxicity"] = sc
            p.methods["toxicity"] = tox_engine
            flagged = sc >= soft_thr
            discard = hard_thr is not None and sc >= hard_thr
        else:
            p.metrics["toxicity"] = None
            p.metrics["toxicity_proxy"] = _tox_score(p.seq)
            p.methods["toxicity"] = "ATLANDI (ToxinPred3 yok; vekil yalnız anotasyon)"
            flagged = discard = False
        p.metrics["toxicity_flag"] = flagged
        if flagged and not discard:
            # yumuşak ceza — scoring.score() bu katsayıyı uygular
            p.metrics["toxicity_penalty_factor"] = penalty
            n_toxic_flag += 1
        if discard:
            p.passed = False
            p.notes.append(f"toksik — yüksek güven (skor ≥ {hard_thr}) → elendi")
            n_toxic_discard += 1
        else:
            survivors.append(p)

    summary = {
        "girdi_peptit": n_start,
        "alerjenite_sonrasi": len(after_allergen), "alerjen_elenen": n_allergen,
        "toksisite_sonrasi": len(survivors),
        "toksik_isaretlenen": n_toxic_flag, "toksik_elenen": n_toxic_discard,
        "alerjenite": "GERÇEK (FAO/WHO 6-mer)" if alg_avail else "ATLANDI",
        "toksisite": tox_engine if tox_avail else "ATLANDI",
        "toksisite_modu": (f"yumuşak işaretleme ≥{soft_thr}"
                           + (f", eleme ≥{hard_thr}" if hard_thr is not None else "")),
        "toksisite_kalibrasyon": "benchmark/results/filter_calibration.json",
    }
    return survivors, summary
