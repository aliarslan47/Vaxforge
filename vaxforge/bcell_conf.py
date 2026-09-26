"""Yapı-tabanlı konformasyonel/yüzey B-hücre epitop öngörüsü.

BepiPred lineer B-epitopu verir; gerçek B-epitoplarının çoğu 3B'de yüzeyde
kümelenen SÜREKSIZ kalıntılardır. Burada ESMFold yapısından RSA (yüzey
erişilebilirliği) + Parker hidrofilisitesi birleştirilerek yüzeye açık, hidrofilik
kalıntı yamaları (patch) konformasyonel B-epitop adayı olarak işaretlenir.

DÜRÜSTLÜK: bu DiscoTope-BENZERİ bir öngörüdür — resmi DiscoTope-3 modeli DEĞİL,
deneysel değil. RSA + hidrofilisite kombinasyonuna dayanan yapısal-yüzey heuristiği.
"""

from __future__ import annotations

# Parker hidrofilisite ölçeği (yüksek = hidrofilik = yüzey/epitop eğilimli).
_PARKER = {
    "A": 2.1, "C": 1.4, "D": 10.0, "E": 7.8, "F": -9.2, "G": 5.7, "H": 2.1,
    "I": -8.0, "K": 5.7, "L": -9.2, "M": -4.2, "N": 7.0, "P": 2.1, "Q": 6.0,
    "R": 4.2, "S": 6.5, "T": 5.2, "V": -3.7, "W": -10.0, "Y": -1.9,
}
MIN_PATCH = 5     # bir yüzey B-yamasında en az kalıntı


def _norm_parker(aa: str) -> float:
    """Parker skorunu 0-1'e getir (min -10, max 10)."""
    v = _PARKER.get(aa.upper(), 0.0)
    return (v + 10.0) / 20.0


def predict(seq: str, rsa_list: list[dict]) -> dict | None:
    """(dizi, RSA listesi) → konformasyonel yüzey B-epitop yamaları.

    Kalıntı skoru = 0.6·RSA + 0.4·Parker(norm). Yüzeyde (exposed) VE skor≥0.5
    olan ardışık ≥MIN_PATCH kalıntı → bir yama. Hata/uyumsuz → None (dürüst).
    """
    if not seq or not rsa_list:
        return None
    n = min(len(seq), len(rsa_list))
    if n < MIN_PATCH:
        return None
    scores = []
    flags = []
    for i in range(n):
        rsa = float(rsa_list[i].get("rsa", 0.0))
        exposed = bool(rsa_list[i].get("exposed"))
        s = 0.6 * min(1.0, rsa) + 0.4 * _norm_parker(seq[i])
        scores.append(round(s, 3))
        flags.append(exposed and s >= 0.5)
    # ardışık bayraklı kalıntıları yamalara topla
    patches = []
    i = 0
    while i < n:
        if flags[i]:
            j = i
            while j < n and flags[j]:
                j += 1
            if j - i >= MIN_PATCH:
                patches.append({
                    "start": i + 1, "end": j,
                    "length": j - i,
                    "peptide": seq[i:j],
                    "mean_score": round(sum(scores[i:j]) / (j - i), 3),
                })
            i = j
        else:
            i += 1
    return {
        "n_patches": len(patches),
        "patches": patches,
        "residues_in_patches": sum(p["length"] for p in patches),
        "method": "yapı-tabanlı yüzey B (DiscoTope-benzeri; RSA+Parker) — öngörü, resmi DiscoTope değil",
    }
