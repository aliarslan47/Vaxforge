"""Adım 7 — Çok-epitoplu mRNA aşı konstrüktü (tam protokol).

Protein katmanı: adjuvan + (EAAAK) + epitoplar (CTL:AAY, HTL:GPGPG, B:KK).
mRNA katmanı: insan kodon optimizasyonu + 5'UTR/Kozak + start + stop + 3'UTR + polyA.
Doğrulama: GC%, CAI (insan), ProtParam (kararsızlık/pI/MW). RNAfold MFE gerçek
ViennaRNA gerektirir; yoksa 'çalıştırılmadı' işaretlenir.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .models import Peptide
from .sequtils import protparam, sanitize

# İnsan için tercih edilen kodon + göreli adaptasyon ağırlığı (w) — kompakt tablo
HUMAN_CODON = {
    "A": ("GCC", 1.0), "R": ("CGG", 1.0), "N": ("AAC", 1.0), "D": ("GAC", 1.0),
    "C": ("TGC", 1.0), "Q": ("CAG", 1.0), "E": ("GAG", 1.0), "G": ("GGC", 1.0),
    "H": ("CAC", 1.0), "I": ("ATC", 1.0), "L": ("CTG", 1.0), "K": ("AAG", 1.0),
    "M": ("ATG", 1.0), "F": ("TTC", 1.0), "P": ("CCC", 1.0), "S": ("AGC", 1.0),
    "T": ("ACC", 1.0), "W": ("TGG", 1.0), "Y": ("TAC", 1.0), "V": ("GTG", 1.0),
}
# insan β-defensin-2 (adjuvan) kısmı — temsili
ADJUVANT = "GIGDPVTCLKSGAICHPVFCPRRYKQIGTCGLPGTKCCKKP"
FIVE_UTR = "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC"  # +Kozak(GCCACC)
THREE_UTR = "GCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCC"
POLY_A = "A" * 30


@dataclass
class Construct:
    protein: str
    mrna: str
    parts: list[str] = field(default_factory=list)
    metrics: dict = field(default_factory=dict)
    methods: dict = field(default_factory=dict)


def _optimize(protein: str) -> tuple[str, float, float]:
    """İnsan kodonlarına çevir; (dna, gc%, CAI) döndür."""
    s = sanitize(protein)
    codons = [HUMAN_CODON[c][0] for c in s if c in HUMAN_CODON]
    dna = "".join(codons)
    gc = 100 * sum(1 for b in dna if b in "GC") / len(dna) if dna else 0
    ws = [HUMAN_CODON[c][1] for c in s if c in HUMAN_CODON]
    cai = round(pow(2.718, sum(__import__("math").log(w) for w in ws) / len(ws)), 3) if ws else 0
    return dna, round(gc, 1), cai


def build(top_peptides: list[Peptide], n_ctl=4, n_htl=4, n_bcell=6) -> Construct:
    """En iyi adaylardan sınıf-başına seçip dengeli çok-epitoplu konstrükt kurar.

    top_peptides adaylık puanına göre sıralı gelir; her sınıftan en iyileri alırız.
    """
    passed = [p for p in top_peptides if p.passed]
    ctl = [p.seq for p in passed if p.kind == "MHC-I"][:n_ctl]
    htl = [p.seq for p in passed if p.kind == "MHC-II"][:n_htl]
    bcell = [p.seq for p in passed if p.kind == "B"][:n_bcell]

    parts = [f"adjuvan(β-defensin)={ADJUVANT}"]
    chunks = [ADJUVANT, "EAAAK"]
    if ctl:
        chunks.append("AAY".join(ctl)); parts.append(f"CTL(AAY)×{len(ctl)}")
    if htl:
        chunks.append("GPGPG" + "GPGPG".join(htl)); parts.append(f"HTL(GPGPG)×{len(htl)}")
    if bcell:
        chunks.append("KK" + "KK".join(bcell)); parts.append(f"B(KK)×{len(bcell)}")
    protein = "".join(chunks)

    cds, gc, cai = _optimize(protein)
    mrna = FIVE_UTR + "ATG" + cds + "TAA" + THREE_UTR + POLY_A

    pp = protparam(protein)
    c = Construct(protein=protein, mrna=mrna, parts=parts)
    c.metrics = {
        "protein_len": len(protein), "mrna_len": len(mrna),
        "gc_percent": gc, "cai_human": cai,
        "instability": pp.get("instability"), "pI": pp.get("pI"),
        "mw_kda": round(pp.get("mw", 0) / 1000, 1),
        "n_ctl": len(ctl), "n_htl": len(htl), "n_bcell": len(bcell),
        "rnafold_mfe": None,
    }
    c.methods = {
        "codon_opt": "insan tercih kodonları (kompakt tablo)",
        "cai": "insan w-ağırlıkları",
        "rnafold_mfe": "ÇALIŞTIRILMADI (ViennaRNA yok)",
        "protparam": "Biopython ProtParam (gerçek)",
    }
    return c
