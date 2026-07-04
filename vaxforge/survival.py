"""Adım 5 — Sağ kalım elemesi (filtre kısmı).

Sert filtreler: alerjenite (non-alerjen olmalı), toksisite.
Yumuşak: konservasyon (suş verisi yoksa 'çalıştırılmadı'), popülasyon/tür kapsamı.

Gerçek araçlar: AllerTOP/AllergenFP, ToxinPred, IEDB Population Coverage.
Buradakiler HEURİSTİK vekillerdir ve 'method' ile etiketlenir.
(2. aşama olan docking/MD ağır adımdır; bu makinede 'deferred'.)
"""

from __future__ import annotations

from . import allergen, toxinpred
from .config_loader import ResolvedTool
from .models import Peptide
from .sequtils import KD, sanitize

# kaba alerjen ilişkili / toksik motif ipuçları (yalnız demonstrasyon)
_TOX_HINT = set("CKR")   # sistein-zengin + yüksek pozitif yük kaba toksisite işareti


def _allergen_prob(pep: str) -> float:
    s = sanitize(pep)
    if not s:
        return 1.0
    aromatic = sum(1 for c in s if c in "FWY") / len(s)
    hydro = sum(1 for c in s if KD.get(c, 0) > 1.5) / len(s)
    return round(min(1.0, 0.3 + 0.6 * aromatic + 0.3 * hydro), 3)


def _tox_score(pep: str) -> float:
    s = sanitize(pep)
    if not s:
        return 0.0
    cys = sum(1 for c in s if c == "C") / len(s)
    pos = sum(1 for c in s if c in "KR") / len(s)
    return round(2 * (0.6 * cys + 0.5 * pos) - 0.4, 3)   # ~ SVM skoru vekili


def run(peptides: list[Peptide], tools: dict[str, ResolvedTool]) -> tuple[list[Peptide], dict]:
    require_non_allergen = bool(tools["allergenicity"].params["require_non_allergen"].value)
    max_tox = float(tools["toxicity"].params["max_toxicity_score"].value)

    active = [p for p in peptides if p.passed]
    seqs = [p.seq for p in active]

    # GERÇEK araçları varsa toplu çalıştır
    alg_avail = allergen.available()
    tox_avail = toxinpred.available()
    alg_res = allergen.predict(seqs) if alg_avail else {}
    tox_res = toxinpred.predict(seqs) if tox_avail else {}

    survivors = []
    n_allergen = n_toxic = 0
    for p in peptides:
        if not p.passed:
            continue
        # -- alerjenite
        if alg_avail and p.seq in alg_res:
            is_allergen = alg_res[p.seq]["allergen"]
            p.metrics["allergen"] = is_allergen
            p.metrics["allergen_match"] = alg_res[p.seq]["match"]
            p.methods["allergen"] = "GERÇEK (FAO/WHO 6-mer, UniProt allergen DB)"
        else:
            ap = _allergen_prob(p.seq)
            p.metrics["allergen"] = ap >= 0.6
            p.metrics["allergen_prob"] = ap
            p.methods["allergen"] = "heuristik proxy (AllerTOP yok)"
            is_allergen = ap >= 0.6
        # -- toksisite
        if tox_avail and p.seq in tox_res:
            toxic = tox_res[p.seq]["toxic"]
            p.metrics["toxicity"] = tox_res[p.seq]["ml_score"]
            p.methods["toxicity"] = "GERÇEK (ToxinPred2, RF)"
        else:
            tx = _tox_score(p.seq)
            p.metrics["toxicity"] = tx
            p.methods["toxicity"] = "heuristik proxy (ToxinPred yok)"
            toxic = tx > max_tox

        if require_non_allergen and is_allergen:
            p.passed = False
            p.notes.append("alerjen (elendi)")
            n_allergen += 1
        if toxic:
            p.passed = False
            p.notes.append("toksik (elendi)")
            n_toxic += 1

        p.metrics.setdefault("conservation", None)
        p.methods["conservation"] = "ÇALIŞTIRILMADI (suş verisi yok)"

        if p.passed:
            survivors.append(p)
    summary = {"girdi": len(peptides), "hayatta": len(survivors),
               "alerjen_elenen": n_allergen, "toksik_elenen": n_toxic,
               "alerjenite": "GERÇEK (FAO/WHO)" if alg_avail else "proxy",
               "toksisite": "GERÇEK (ToxinPred2)" if tox_avail else "proxy"}
    return survivors, summary
