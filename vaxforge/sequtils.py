"""Dizi/protein özellik yardımcıları (saf-Python + Biopython).

Buradaki fonksiyonlar harici araç gerektirmeden bu makinede koşan
HEURİSTİK karşılıklardır. Gerçek araçlar (DeepTMHMM, SignalP, VaxiJen...)
takıldığında ilgili modül bunların yerine gerçek sonucu kullanır.
"""

from __future__ import annotations

from Bio.SeqUtils.ProtParam import ProteinAnalysis

_STD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hidropati indeksi
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Parker hidrofilisite (B-hücre epitobu eğilimi için)
PARKER = {
    "A": 2.1, "R": 4.2, "N": 7.0, "D": 10.0, "C": 1.4, "Q": 6.0, "E": 7.8,
    "G": 5.7, "H": 2.1, "I": -8.0, "L": -9.2, "K": 5.7, "M": -4.2, "F": -9.2,
    "P": 2.1, "S": 6.5, "T": 5.2, "W": -10.0, "Y": -1.9, "V": -3.7,
}

# Kolaskar-Tongaonkar (1990) antijenite eğilimi — yayınlanmış yarı-ampirik ölçek.
# Klasik B-hücre epitop antijenite yöntemi (IEDB'de de mevcut). Eşik: ortalama ≥ 1.0
# (veya protein ortalaması) ise antijenik.
KOLASKAR = {
    "A": 1.064, "R": 0.873, "N": 0.776, "D": 0.866, "C": 1.412, "Q": 1.015,
    "E": 0.851, "G": 0.874, "H": 1.105, "I": 1.152, "L": 1.250, "K": 0.930,
    "M": 0.826, "F": 1.091, "P": 1.064, "S": 1.012, "T": 0.909, "W": 0.893,
    "Y": 1.161, "V": 1.383,
}


def kolaskar_mean(pep: str) -> float:
    s = sanitize(pep)
    if not s:
        return 0.0
    return sum(KOLASKAR.get(c, 1.0) for c in s) / len(s)


def parker_mean(pep: str) -> float:
    s = sanitize(pep)
    if not s:
        return 0.0
    return sum(PARKER.get(c, 0.0) for c in s) / len(s)


def sanitize(seq: str) -> str:
    """Standart-dışı harfleri temizler (ProteinAnalysis için)."""
    s = seq.strip().upper().replace("*", "")
    return "".join(c for c in s if c in _STD_AA)


def protparam(seq: str) -> dict:
    """ProtParam özet özellikleri (kararsızlık, pI, GRAVY, aromatiklik...)."""
    s = sanitize(seq)
    if len(s) < 5:
        return {}
    pa = ProteinAnalysis(s)
    return {
        "length": len(s),
        "mw": round(pa.molecular_weight(), 1),
        "pI": round(pa.isoelectric_point(), 2),
        "instability": round(pa.instability_index(), 1),
        "gravy": round(pa.gravy(), 3),
        "aromaticity": round(pa.aromaticity(), 3),
    }


def hydropathy_profile(seq: str, window: int = 19) -> list[float]:
    s = sanitize(seq)
    if len(s) < window:
        return []
    vals = [KD.get(c, 0.0) for c in s]
    out = []
    for i in range(len(s) - window + 1):
        out.append(sum(vals[i:i + window]) / window)
    return out


def count_tm_helices(seq: str, window: int = 19, thresh: float = 1.6) -> int:
    """Kyte-Doolittle pencere ortalaması eşiği aşan bölge sayısı ≈ TM heliks."""
    prof = hydropathy_profile(seq, window)
    count, inside = 0, False
    for v in prof:
        if v >= thresh and not inside:
            count += 1
            inside = True
        elif v < thresh:
            inside = False
    return count


def signal_peptide_prob(seq: str) -> float:
    """N-uç sinyal peptidi HEURİSTİĞİ: ilk ~30 aa'da hidrofobik h-bölge var mı?

    Gerçek: SignalP/DeepSig. Bu proxy 0-1 arası kaba bir olasılık verir.
    """
    s = sanitize(seq)
    if len(s) < 20:
        return 0.0
    nterm = s[:30]
    # h-bölge: 7-15 arası hidrofobik pencere
    best = 0.0
    for i in range(3, min(len(nterm) - 7, 18)):
        win = nterm[i:i + 8]
        h = sum(1 for c in win if KD.get(c, 0) > 1.5) / 8
        best = max(best, h)
    # n-bölge pozitif yük bonusu
    pos = sum(1 for c in nterm[:5] if c in "KR") / 5
    return round(min(1.0, 0.7 * best + 0.3 * pos), 3)


def predict_localization(seq: str) -> str:
    """Kaba lokalizasyon HEURİSTİĞİ (gerçek: DeepLoc/PSORTb).

    Sinyal peptidi + TM sayısına göre yüzey/salgı vs sitoplazma tahmini.
    """
    sp = signal_peptide_prob(seq)
    tm = count_tm_helices(seq)
    if tm >= 2:
        return "membrane"
    if sp >= 0.5 and tm <= 1:
        return "secreted" if tm == 0 else "outer_membrane"
    if tm == 1:
        return "outer_membrane"
    return "cytoplasm"


def antigenicity_proxy(seq: str) -> float:
    """VaxiJen YERİNE fizikokimyasal antijenite proxy'si (0-1).

    Gerçek VaxiJen ACC dönüşümü + SVM kullanır. Bu, hidrofilisite +
    yük + esneklik karışımından türetilen kaba bir vekildir; SADECE
    demonstrasyon amaçlıdır, gerçek VaxiJen skoru değildir.
    """
    s = sanitize(seq)
    if not s:
        return 0.0
    hyd = sum(PARKER.get(c, 0) for c in s) / len(s)      # yüksek = hidrofilik
    charged = sum(1 for c in s if c in "DEKR") / len(s)
    flex = sum(1 for c in s if c in "GSPTAN") / len(s)
    raw = 0.4 + 0.02 * hyd + 0.4 * charged + 0.3 * flex
    return round(max(0.0, min(1.0, raw)), 3)
