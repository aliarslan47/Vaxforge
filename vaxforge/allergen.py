"""GERÇEK alerjenite tahmini — FAO/WHO biyoinformatik ölçütü (çevrimdışı).

Kısa epitoplar için FAO/WHO Codex ölçütü: bir peptit, bilinen bir allerjenle
ardışık ≥6 aynı amino asit paylaşıyorsa potansiyel alerjen sayılır.

Bilinen allerjen seti: UniProt KW-0020 (allergen anahtar kelimesi, ~1021 protein),
tools/db/allergens.fasta.gz. Tüm 6-mer'leri bir kez çıkarıp önbelleğe (pickle)
alır. Set yoksa available()=False -> survival.py proxy'ye düşer.
"""

from __future__ import annotations

import functools
import gzip
import pickle
from pathlib import Path

_DB = Path(__file__).resolve().parent.parent / "tools" / "db" / "allergens.fasta.gz"
_CACHE = Path(__file__).resolve().parent.parent / "tools" / "db" / "allergen_6mers.pkl"
_K = 6


def available() -> bool:
    return _DB.exists()


@functools.lru_cache(maxsize=1)
def _kmers() -> frozenset:
    if _CACHE.exists():
        try:
            return frozenset(pickle.loads(_CACHE.read_bytes()))
        except Exception:
            pass
    kmers: set[str] = set()
    if not _DB.exists():
        return frozenset()
    with gzip.open(_DB, "rt", encoding="latin-1") as fh:
        seq = []
        for line in fh:
            if line.startswith(">"):
                _add_kmers("".join(seq), kmers)
                seq = []
            else:
                seq.append(line.strip())
        _add_kmers("".join(seq), kmers)
    try:
        _CACHE.write_bytes(pickle.dumps(kmers))
    except Exception:
        pass
    return frozenset(kmers)


def _add_kmers(seq: str, out: set) -> None:
    s = seq.upper()
    for i in range(len(s) - _K + 1):
        out.add(s[i:i + _K])


def predict(peptides: list[str]) -> dict[str, dict]:
    """Peptit -> {allergen, match}. 6-mer allerjen eşleşmesi (FAO/WHO)."""
    kmers = _kmers()
    if not kmers:
        return {}
    out: dict[str, dict] = {}
    for pep in set(peptides):
        s = pep.upper()
        hit = ""
        for i in range(len(s) - _K + 1):
            if s[i:i + _K] in kmers:
                hit = s[i:i + _K]
                break
        out[pep] = {"allergen": bool(hit), "match": hit}
    return out
