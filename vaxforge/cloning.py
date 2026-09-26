"""In silico klonlama — kodon optimizasyonu (CAI) + restriksiyon + vektör önerisi.

MEV dizisi rekombinant olarak üretilecekse konak (E. coli) kodon kullanımına
uyarlanmalı: nadir kodonlar ekspresyonu düşürür. Bu modül MEV proteinini seçili
konağın EN SIK kodonlarıyla ters-translasyona çevirir, CAI (Codon Adaptation Index;
Sharp & Li 1987) ve GC% hesaplar, iç restriksiyon bölgelerini tarar (klonlamayı
bozabilir) ve uygun ekspresyon vektörü/etiket önerir.

Kaynak: E. coli K-12 kodon kullanımı (Kazusa). Saf Python + Biopython.Restriction —
dış araç yok. DÜRÜSTLÜK: CAI/GC hesaplanır; 'optimallik' konağa görelidir, deneysel
ekspresyon garantisi değildir.
"""

from __future__ import annotations

# E. coli K-12 kodon kullanım fraksiyonu (aminoasit-içi göreli), Kazusa.
# {aa: {codon: fraction}}. CAI ağırlığı w = fraction / max(fraction ailesi).
_ECOLI = {
    "A": {"GCG": 0.36, "GCC": 0.27, "GCA": 0.21, "GCT": 0.16},
    "R": {"CGC": 0.40, "CGT": 0.38, "CGG": 0.10, "CGA": 0.06, "AGA": 0.04, "AGG": 0.02},
    "N": {"AAC": 0.55, "AAT": 0.45},
    "D": {"GAT": 0.63, "GAC": 0.37},
    "C": {"TGC": 0.56, "TGT": 0.44},
    "Q": {"CAG": 0.65, "CAA": 0.35},
    "E": {"GAA": 0.68, "GAG": 0.32},
    "G": {"GGC": 0.40, "GGT": 0.34, "GGG": 0.15, "GGA": 0.11},
    "H": {"CAT": 0.57, "CAC": 0.43},
    "I": {"ATT": 0.51, "ATC": 0.42, "ATA": 0.07},
    "L": {"CTG": 0.50, "TTA": 0.13, "TTG": 0.13, "CTC": 0.10, "CTT": 0.10, "CTA": 0.04},
    "K": {"AAA": 0.76, "AAG": 0.24},
    "M": {"ATG": 1.00},
    "F": {"TTT": 0.57, "TTC": 0.43},
    "P": {"CCG": 0.52, "CCA": 0.19, "CCT": 0.16, "CCC": 0.12},
    "S": {"AGC": 0.28, "TCT": 0.15, "TCC": 0.15, "AGT": 0.15, "TCG": 0.15, "TCA": 0.12},
    "T": {"ACC": 0.44, "ACG": 0.27, "ACT": 0.17, "ACA": 0.13},
    "W": {"TGG": 1.00},
    "Y": {"TAT": 0.57, "TAC": 0.43},
    "V": {"GTG": 0.37, "GTT": 0.26, "GTC": 0.22, "GTA": 0.15},
}
_STOP = "TAA"   # E. coli en sık stop kodonu

# En sık kodon (optimize ters-translasyon için) + CAI ağırlıkları.
_BEST = {aa: max(cod, key=cod.get) for aa, cod in _ECOLI.items()}
_W = {codon: frac / max(fam.values())
      for aa, fam in _ECOLI.items() for codon, frac in fam.items()}

# Klonlamada sık kullanılan enzimler (iç bölge = klonlama çakışması riski).
_CLONING_ENZYMES = ["NdeI", "XhoI", "EcoRI", "BamHI", "HindIII", "NcoI", "SalI", "NotI"]


def available() -> bool:
    """Saf Python + Biopython.Restriction — her zaman kullanılabilir."""
    try:
        import Bio.Restriction  # noqa: F401
        return True
    except Exception:
        return False


def back_translate(protein: str, host: str = "ecoli") -> str:
    """Protein → konağın en sık kodonlarıyla optimize DNA (+ stop)."""
    dna = "".join(_BEST.get(aa.upper(), "NNN") for aa in protein if aa.upper() in _BEST)
    return dna + _STOP


def cai(dna: str) -> float:
    """Codon Adaptation Index (Sharp & Li 1987) — w'lerin geometrik ortalaması."""
    import math
    codons = [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
    ws = [_W[c] for c in codons if c in _W]
    if not ws:
        return 0.0
    return round(math.exp(sum(math.log(w) for w in ws) / len(ws)), 3)


def _gc(dna: str) -> float:
    if not dna:
        return 0.0
    return round(100.0 * sum(1 for c in dna if c in "GC") / len(dna), 1)


def _restriction_scan(dna: str) -> dict[str, list[int]]:
    """İç restriksiyon bölgeleri (klonlama çakışması). {enzim: [pozisyon...]}."""
    try:
        from Bio.Restriction import RestrictionBatch
        from Bio.Seq import Seq
        rb = RestrictionBatch(_CLONING_ENZYMES)
        hits = rb.search(Seq(dna))
        return {str(enz): [int(p) for p in sites] for enz, sites in hits.items() if sites}
    except Exception:
        return {}


def run(protein_seq: str, host: str = "ecoli") -> dict | None:
    """MEV proteini → optimize DNA + CAI + GC + restriksiyon + vektör önerisi.

    Döner: {host, protein_len, dna, dna_len_bp, gc_percent, cai, restriction_sites,
            n_internal_sites, suggested_vector, tag, note}. Hata → None.
    """
    if not protein_seq or not available():
        return None
    try:
        dna = back_translate(protein_seq, host)
        sites = _restriction_scan(dna)
        n_internal = sum(len(v) for v in sites.values())
        # NdeI/XhoI iç bölgesi yoksa pET-28a(+) His-tag klonlaması temiz.
        clean_ndex = "NdeI" not in sites and "XhoI" not in sites
        vector = "pET-28a(+)"
        tag = "N-terminal His6 (NdeI–XhoI)" if clean_ndex else \
              "pET-28a(+) — NdeI/XhoI iç bölgesi VAR, alternatif kesim/gen-sentezi gerekir"
        return {
            "available": True,
            "host": "Escherichia coli K-12",
            "method": "kodon-optimizasyon (Kazusa E. coli en sık kodon) + CAI (Sharp&Li 1987)",
            "protein_len": len([a for a in protein_seq if a.upper() in _BEST]),
            "dna": dna,
            "dna_len_bp": len(dna),
            "gc_percent": _gc(dna),
            "cai": cai(dna),
            "restriction_sites": sites,
            "n_internal_sites": n_internal,
            "suggested_vector": vector,
            "tag": tag,
            "note": ("CAI>0.8 iyi ekspresyon beklentisi (E. coli); GC %30-70 arası tercih. "
                     "İç restriksiyon bölgeleri klasik kesim-ligasyonu bozar → gen sentezi/"
                     "alternatif enzim. Bu in silico öngörüdür, deneysel ekspresyon garantisi değil."),
        }
    except Exception:
        return None
