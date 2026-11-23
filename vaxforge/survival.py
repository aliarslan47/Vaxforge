"""Adım 5 — Sağ kalım elemesi (filtre kısmı).

Sert filtreler: alerjenite (non-alerjen olmalı), toksisite.
Yumuşak: konservasyon (suş verisi yoksa 'çalıştırılmadı'), popülasyon/tür kapsamı.

Gerçek araçlar: AllerTOP/AllergenFP, ToxinPred, IEDB Population Coverage.
Buradakiler HEURİSTİK vekillerdir ve 'method' ile etiketlenir.
(Yapısal kararlılık/docking-MD ikinci aşaması şimdilik çıkarıldı.)
"""

from __future__ import annotations

from . import allergen, toxinpred
from .config_loader import ResolvedTool
from .i18n import t
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
    tox_res = toxinpred.predict(seqs, threshold=max_tox) if tox_avail else {}

    n_start = len(active)

    # --- Aşama 1: ALERJENİTE elemesi ---------------------------------------
    after_allergen, n_allergen = [], 0
    for p in active:
        if alg_avail and p.seq in alg_res:
            is_allergen = alg_res[p.seq]["allergen"]
            p.metrics["allergen"] = is_allergen
            p.metrics["allergen_match"] = alg_res[p.seq]["match"]
            p.methods["allergen"] = "GERÇEK (FAO/WHO 6-mer, UniProt allergen DB)"
        else:
            ap = _allergen_prob(p.seq)
            is_allergen = ap >= 0.6
            p.metrics["allergen"] = is_allergen
            p.metrics["allergen_prob"] = ap
            p.methods["allergen"] = "heuristik proxy (AllerTOP yok)"
        if require_non_allergen and is_allergen:
            p.passed = False
            p.notes.append("alerjen (elendi)")
            n_allergen += 1
        else:
            after_allergen.append(p)

    # --- Aşama 2: TOKSİSİTE elemesi ----------------------------------------
    #
    # UZUNLUK ALANI UYARISI: ToxinPred2'nin altı referans setinin tamamında en
    # kısa dizi 35 aa. Bizim huniye verdiğimiz MHC epitopları 8-15 aa, yani
    # aracın DOĞRULANMADIĞI uzunluk alanı. Bu ortamda ölçüldü (Model 1, eşik
    # 0.6, referans setlerden rastgele 9-mer'ler, n=100+100):
    #     tam dizi (>=35 aa) AUC 1.00  /  9-mer AUC 0.69
    #     eşik 0.6'da toksik-OLMAYAN proteinden gelen 9-mer'lerin 66/100'ü
    #     "toksik" ilan ediliyor.
    # Hiçbir eşik iyi denge vermiyor (0.7 → seçicilik .77/duyarlılık .49).
    #
    # Bu yüzden alan DIŞI peptitte skor HESAPLANIR ve raporlanır ama SERT
    # FİLTRE uygulanmaz. Aksi hâlde güvenli epitopların üçte ikisi sessizce
    # elenirdi. (Araç daha önce bozuk olduğu için bu kalibrasyon sorunu
    # görünmüyordu; bkz. DURUM.md bulgu 22.)
    survivors, n_toxic, n_tox_ood = [], 0, 0
    for p in after_allergen:
        in_domain = toxinpred.in_validated_domain(p.seq)
        if tox_avail and p.seq in tox_res:
            toxic = tox_res[p.seq]["toxic"]
            p.metrics["toxicity"] = tox_res[p.seq]["ml_score"]
            p.methods["toxicity"] = "GERÇEK (ToxinPred2, RF)"
            if not in_domain:
                # ölçülmüş ama filtreye SOKULMAZ — dürüst işaret
                toxic = False
                n_tox_ood += 1
                p.methods["toxicity"] = (
                    f"GERÇEK (ToxinPred2, RF) — UYARI: peptit "
                    f"{len(p.seq)} aa < {toxinpred.MIN_VALIDATED_LENGTH} aa "
                    "doğrulanmış alan; skor bilgi amaçlı, eleme yapılmadı")
                p.notes.append(t("tr", "tox_out_of_domain"))
                p.notes.append(t("en", "tox_out_of_domain"))
        else:
            tx = _tox_score(p.seq)
            toxic = tx > max_tox
            p.metrics["toxicity"] = tx
            p.methods["toxicity"] = "heuristik proxy (ToxinPred yok)"
        if toxic:
            p.passed = False
            p.notes.append("toksik (elendi)")
            n_toxic += 1
        else:
            survivors.append(p)

    summary = {
        "girdi_peptit": n_start,
        "alerjenite_sonrasi": len(after_allergen), "alerjen_elenen": n_allergen,
        "toksisite_sonrasi": len(survivors), "toksik_elenen": n_toxic,
        "toksisite_alan_disi": n_tox_ood,
        "alerjenite": "GERÇEK (FAO/WHO)" if alg_avail else "proxy",
        "toksisite": "GERÇEK (ToxinPred2)" if tox_avail else "proxy",
    }
    return survivors, summary
