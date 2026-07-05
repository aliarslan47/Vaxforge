"""Adım 4 (MHC-I sonrası) — CD8+ T-hücre immünojenite skoru.

GERÇEK model: IEDB Class-I Immunogenicity (Calis ve ark. 2013, PLoS Comput Biol).
Peptidin immünojenitesini amino asit kullanımının pozisyona-bağlı zenginleşmesiyle
tahmin eder. MHC bağlanmasından SONRA uygulanır: bağlanan bir peptidin gerçekten
T-hücre yanıtı doğurma olasılığını ölçer (bağlanma ≠ immünojenite).

Sabitler (immunoscale, immunoweight, allel maskeleri) IEDB'nin resmi aracından
(tools/immunogenicity/predict_immunogenicity.py, IEDB_Immunogenicity-3.0) BİREBİR
alınmıştır — proxy değil, atıflanabilir referans model. Skor > 0 → immünojenik
eğilim; < 0 → tersi. `score_norm` sıralama için 0-1'e eşlenmiş monoton hâli.
"""

from __future__ import annotations

import math

# AA immünojenite ölçeği (Calis 2013)
_IMMUNOSCALE = {
    "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110,
    "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036, "M": -0.570, "N": -0.021,
    "P": -0.036, "Q": -0.376, "R": 0.168, "S": -0.537, "T": 0.126, "V": 0.134,
    "W": 0.719, "Y": -0.012,
}
# 9-mer pozisyon ağırlıkları (Calis 2013)
_IMMUNOWEIGHT = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

# Allele -> varsayılan anchor maske pozisyonları (1-tabanlı) — IEDB allele_dict
_ALLELE_MASK = {
    "H-2-Db": "2,5,9", "H-2-Dd": "2,3,5", "H-2-Kb": "2,3,9", "H-2-Kd": "2,5,9",
    "H-2-Kk": "2,8,9", "H-2-Ld": "2,5,9", "HLA-A0101": "2,3,9", "HLA-A0201": "1,2,9",
    "HLA-A0202": "1,2,9", "HLA-A0203": "1,2,9", "HLA-A0206": "1,2,9", "HLA-A0211": "1,2,9",
    "HLA-A0301": "1,2,9", "HLA-A1101": "1,2,9", "HLA-A2301": "2,7,9", "HLA-A2402": "2,7,9",
    "HLA-A2601": "1,2,9", "HLA-A2902": "2,7,9", "HLA-A3001": "1,3,9", "HLA-A3002": "2,7,9",
    "HLA-A3101": "1,2,9", "HLA-A3201": "1,2,9", "HLA-A3301": "1,2,9", "HLA-A6801": "1,2,9",
    "HLA-A6802": "1,2,9", "HLA-A6901": "1,2,9", "HLA-B0702": "1,2,9", "HLA-B0801": "2,5,9",
    "HLA-B1501": "1,2,9", "HLA-B1502": "1,2,9", "HLA-B1801": "1,2,9", "HLA-B2705": "2,3,9",
    "HLA-B3501": "1,2,9", "HLA-B3901": "1,2,9", "HLA-B4001": "1,2,9", "HLA-B4002": "1,2,9",
    "HLA-B4402": "2,3,9", "HLA-B4403": "2,3,9", "HLA-B4501": "1,2,9", "HLA-B4601": "1,2,9",
    "HLA-B5101": "1,2,9", "HLA-B5301": "1,2,9", "HLA-B5401": "1,2,9", "HLA-B5701": "1,2,9",
    "HLA-B5801": "1,2,9",
}

METHOD = "GERÇEK (IEDB Class-I Immunogenicity, Calis 2013)"


def available() -> bool:
    return True  # saf-Python, harici bağımlılık yok


def _allele_key(allele: str | None) -> str | None:
    """NetMHCpan allel adını IEDB anahtar biçimine indirger (HLA-A*02:01 -> HLA-A0201)."""
    if not allele:
        return None
    return allele.replace("*", "").replace(":", "").strip()


def score(peptide: str, allele: str | None = None) -> float | None:
    """Calis 2013 immünojenite ham skoru. Geçersiz peptit -> None.

    Maskeleme: allel bilinip IEDB tablosunda varsa ona özel anchor pozisyonları;
    yoksa varsayılan (1, 2, C-terminal) maskelenir (aracın default davranışı).
    """
    pep = (peptide or "").upper()
    peplen = len(pep)
    if peplen < 3 or any(a not in _IMMUNOSCALE for a in pep):
        return None

    cterm = peplen - 1
    key = _allele_key(allele)
    mask_str = _ALLELE_MASK.get(key) if key else None
    if mask_str:
        mask_num = [int(x) - 1 for x in mask_str.split(",")]
    else:
        mask_num = [0, 1, cterm]   # varsayılan: 1, 2, C-terminal

    if peplen > 9:
        pepweight = _IMMUNOWEIGHT[:5] + ((peplen - 9) * [0.30]) + _IMMUNOWEIGHT[5:]
    else:
        pepweight = _IMMUNOWEIGHT

    s = 0.0
    for count, aa in enumerate(pep):
        if count in mask_num or count >= len(pepweight):
            continue
        s += pepweight[count] * _IMMUNOSCALE[aa]
    return round(s, 5)


def score_norm(raw: float | None) -> float:
    """Ham skoru sıralama için 0-1'e eşler (monoton sigmoid; 0 -> 0.5)."""
    if raw is None:
        return 0.0
    return round(1.0 / (1.0 + math.exp(-3.0 * raw)), 4)
