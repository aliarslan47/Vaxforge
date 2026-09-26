"""Glikozilasyon taraması — glikan-maskeli B-hücre epitop filtresi.

Patojen yüzey proteinleri N-bağlı glikozilasyon bölgelerine (N-X-S/T, X≠Pro)
taktıkları glikanlarla epitopları antikordan MASKELER (glikan kalkanı; ör. HIV
Env, SARS-CoV-2 spike). Bir B-hücre epitobu bir sequon ile örtüşüyorsa in vivo
glikanla kaplı olabilir → antikor erişemez → zayıf aşı hedefi.

DÜRÜSTLÜK: N-X-S/T dizi motifi KESİN kuraldır (sequon), ama gerçek glikan
işgali deneysel doğrulanmalı; O-glikozilasyon burada yalnız Ser/Thr-yoğun
HEURİSTİK hotspot'tur (kesin değil). Bu bir ÖNGÖRÜDÜR — bölge yumuşak işaretlenir,
elenmez (aday kaybı yok). Bkz. yapısal RSA (structure.py) ile birlikte 'erişilebilirlik'.
"""

from __future__ import annotations

import re

# N-bağlı sequon: N, ardından Pro olmayan herhangi bir kalıntı, ardından S veya T.
_N_SEQUON = re.compile(r"N[^P][ST]")
# O-glikozilasyon HEURİSTİĞİ: 5-pencerede ≥3 Ser/Thr → yoğun hotspot (kesin değil).
_OT_WINDOW = 5
_OT_MIN = 3


def available() -> bool:
    """Saf dizi motifi taraması — dış araç gerekmez (her zaman kullanılabilir)."""
    return True


def scan_sequons(seq: str) -> dict:
    """Bir dizide N-glikozilasyon sequon'ları + O-gliko hotspot'ları.

    Döner: {n_sequons: [(pos, motif)...], n_count, o_hotspots: [pos...], o_count}
    (pos 1-tabanlı, motifin N'sinin konumu.)
    """
    s = (seq or "").upper()
    n_sites = [(m.start() + 1, m.group()) for m in _N_SEQUON.finditer(s)]
    o_sites = []
    for i in range(len(s) - _OT_WINDOW + 1):
        win = s[i:i + _OT_WINDOW]
        if sum(1 for c in win if c in "ST") >= _OT_MIN:
            o_sites.append(i + 1)
    return {
        "n_sequons": n_sites, "n_count": len(n_sites),
        "o_hotspots": o_sites, "o_count": len(o_sites),
    }


def annotate_peptides(peptides: list) -> dict:
    """Aday peptitleri glikozilasyon açısından işaretler (YUMUŞAK — elemez).

    Her peptit: metrics['glyco_n_sequons'], metrics['glyco_masked'] (N-sequon içeriyorsa),
    methods['glyco']. B-hücre epitoplarında maskeleme antikor erişimi açısından KRİTİK;
    T-hücre (MHC) epitoplarında bilgi amaçlı işaretlenir. Döner: özet sözlük.
    """
    n_flag_b = 0
    n_bcell = 0
    n_flag_all = 0
    for p in peptides:
        sc = scan_sequons(p.seq)
        masked = sc["n_count"] > 0
        p.metrics["glyco_n_sequons"] = sc["n_count"]
        p.metrics["glyco_o_hotspots"] = sc["o_count"]
        p.metrics["glyco_masked"] = masked
        p.methods["glyco"] = "N-X-S/T sequon taraması (N-gliko); O = Ser/Thr-yoğun heuristik"
        if masked:
            n_flag_all += 1
        if getattr(p, "kind", "") == "B":
            n_bcell += 1
            if masked:
                n_flag_b += 1
    return {
        "available": True,
        "method": "N-X-S/T sequon (kesin motif) + O Ser/Thr-yoğun heuristik — öngörü, deneysel değil",
        "n_candidates": len(peptides),
        "n_masked": n_flag_all,
        "n_bcell": n_bcell,
        "n_bcell_masked": n_flag_b,
        "note": ("Glikan-maskeli B-epitopları in vivo antikora erişilemeyebilir "
                 "(yumuşak uyarı; aday elenmedi). Yapısal RSA ile birlikte değerlendirin."),
    }
