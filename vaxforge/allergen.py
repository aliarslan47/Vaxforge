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
#: kmer -> kaynak alerjen (ilk görülen) — konstrükt uyarısının "hangi alerjene
#: benziyor" sorusunu yanıtlaması için. 142.320 6-mer, ~15 MB pickle.
_SRC_CACHE = Path(__file__).resolve().parent.parent / "tools" / "db" / "allergen_6mer_src.pkl"
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


# ── Konstrükt düzeyi ayrıntılı eşleşme ─────────────────────────────────────
# NEDEN (2026-11-22): MEV konstrüktünün KENDİSİ FAO/WHO 6-mer ölçütüne göre
# alerjen çıkıyordu ama rapor bunu hiçbir uyarıya çevirmiyordu. Adaylarda
# `require_non_allergen` SERT filtre var, yani her epitop temiz; buna karşın
# birleştirilmiş dizi alerjen çıkıyorsa eşleşme ya LİNKER'dan ya da epitop-linker
# BİRLEŞME YERİNDE doğan yeni 6-mer'den geliyor olmalıdır. Bu ayrım eyleme
# dönüşür: birleşme yerinden doğuyorsa linker değiştirilebilir, epitoptan
# geliyorsa epitop seçimi yanlıştır. `predict` yalnız ilk eşleşmeyi ve bir bool
# döndürdüğü için bu ayrım yapılamıyordu.

@functools.lru_cache(maxsize=1)
def _kmer_sources() -> dict:
    """kmer -> kaynak alerjen etiketi (ilk görülen kayıt). DB yoksa boş."""
    if _SRC_CACHE.exists():
        try:
            return pickle.loads(_SRC_CACHE.read_bytes())
        except Exception:
            pass
    src: dict[str, str] = {}
    if not _DB.exists():
        return {}
    def _add(seq: str, label: str) -> None:
        s = seq.upper()
        for i in range(len(s) - _K + 1):
            src.setdefault(s[i:i + _K], label)
    with gzip.open(_DB, "rt", encoding="latin-1") as fh:
        seq, label = [], ""
        for line in fh:
            if line.startswith(">"):
                if seq:
                    _add("".join(seq), label)
                seq = []
                h = line[1:].strip()
                parts = h.split("|")
                acc = parts[1] if len(parts) > 2 else h.split()[0]
                name = h.split(" ", 1)[1] if " " in h else ""
                os_ = name.split("OS=")[1].split("OX=")[0].strip() if "OS=" in name else ""
                label = f"{acc} ({os_})" if os_ else acc
            else:
                seq.append(line.strip())
        if seq:
            _add("".join(seq), label)
    try:
        _SRC_CACHE.write_bytes(pickle.dumps(src))
    except Exception:
        pass
    return src


def matches(seq: str) -> list[dict]:
    """Dizideki TÜM alerjen 6-mer eşleşmeleri: [{kmer, start, source}].

    `predict` ilk eşleşmede durur (aday elemesi için yeterli); konstrükt
    uyarısı ise hepsini ve YERLERİNİ ister.
    """
    kmers = _kmers()
    if not kmers:
        return []
    srcs = _kmer_sources()
    s = (seq or "").upper()
    out = []
    for i in range(len(s) - _K + 1):
        km = s[i:i + _K]
        if km in kmers:
            out.append({"kmer": km, "start": i, "source": srcs.get(km, "")})
    return out
