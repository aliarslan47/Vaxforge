"""Yüzey erişilebilirliği (RSA) — PDB'den Shrake-Rupley (Biopython).

Kalıntı-başı çözücü-erişilebilir yüzey alanı (SASA) hesaplanır, kalıntı-tipine
göre teorik maksimuma (Tien 2013) bölünerek RELATİF erişilebilirlik (RSA, 0-1)
elde edilir. RSA>0.20 → yüzeyde (antikora erişilebilir); gömülü kalıntılar epitop
olarak zayıftır. Ters aşıda: gömülü tahmini B-epitopları eleyip yüzeydekileri öne çıkarır.
"""

from __future__ import annotations

from io import StringIO

# Tien et al. 2013 (theoretical) — kalıntı maksimum ASA (Å²).
_MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
SURFACE_THRESHOLD = 0.20


def available() -> bool:
    try:
        from Bio.PDB.SASA import ShrakeRupley  # noqa: F401
        return True
    except Exception:
        return False


def rsa_from_pdb(pdb_text: str) -> list[dict] | None:
    """PDB metni → kalıntı-başı [{resnum, resname, rsa, exposed}]. Hata → None."""
    if not available() or not pdb_text:
        return None
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("m", StringIO(pdb_text))
        model = next(iter(structure))
        ShrakeRupley().compute(model, level="R")
        out: list[dict] = []
        for chain in model:
            for res in chain:
                rn = res.resname.strip().upper()
                if rn not in _MAX_ASA:
                    continue
                sasa = float(getattr(res, "sasa", 0.0) or 0.0)
                rsa = sasa / _MAX_ASA[rn]
                out.append({
                    "resnum": res.id[1],
                    "resname": rn,
                    "rsa": round(rsa, 3),
                    "exposed": rsa >= SURFACE_THRESHOLD,
                })
        return out or None
    except Exception:
        return None


def surface_summary(rsa_list: list[dict]) -> dict:
    """RSA listesinden özet: yüzey fraksiyonu + ortalama RSA."""
    if not rsa_list:
        return {}
    n = len(rsa_list)
    exposed = sum(1 for r in rsa_list if r["exposed"])
    return {
        "n_res": n,
        "n_exposed": exposed,
        "percent_surface": round(100.0 * exposed / n, 1),
        "mean_rsa": round(sum(r["rsa"] for r in rsa_list) / n, 3),
    }
