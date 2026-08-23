"""Epitop-seviyesi moleküler mimikri — konak self-peptidine tam eşleşme.

Bir aday epitop, KONAĞIN kendi proteinindeki bir parçayla birebir aynıysa
'moleküler mimikri' = otoimmünite / merkezi-tolerans (T-hücre negatif seçilim)
riski taşır (Vaccine Design Ch4: otoimmün çapraz-reaktif epitopları çıkar).

Protein-seviyesi insan/konak homolojisi funnel'da ZATEN sert filtre; bu modül
GRANULAR (epitop-seviyesi) uyarı katmanıdır: YUMUŞAK — epitopu ELEMEZ, yalnız
işaretler + skorda hafif ceza (kullanıcı kararı). k=9 (MHC-I çekirdek uzunluğu;
paylaşılan 9-mer bir CD8+ epitobu olabilir), config'ten ayarlanır.

Şablon: allergen.py (FAO/WHO 6-mer) — aynı önbellekli k-mer küme deseni. Konak
self-proteomları tools/db/{host}_sprot.fasta.gz (human/bovine/mouse/pig/chicken).
Set yoksa o konak için available()=False → dürüstçe atlanır (uydurma yok).
"""

from __future__ import annotations

import functools
import gzip
import pickle
from pathlib import Path

_DB_DIR = Path(__file__).resolve().parent.parent / "tools" / "db"

# ImmForge tür adı / VaxForge host adı → self-proteom dosya öneki
_HOST_DB = {
    "human": "human", "mouse": "mouse", "bovine": "bovine",
    "pig": "pig", "chicken": "chicken",
}


def db_path(host: str) -> Path:
    return _DB_DIR / f"{_HOST_DB.get(host, host)}_sprot.fasta.gz"


def available(host: str) -> bool:
    return host in _HOST_DB and db_path(host).exists()


@functools.lru_cache(maxsize=8)
def _kmers(host: str, k: int) -> frozenset:
    """Konak self-proteomunun tüm k-mer'leri (önbellekli pickle)."""
    cache = _DB_DIR / f"{_HOST_DB.get(host, host)}_{k}mers.pkl"
    if cache.exists():
        try:
            return frozenset(pickle.loads(cache.read_bytes()))
        except Exception:
            pass
    db = db_path(host)
    if not db.exists():
        return frozenset()
    kmers: set[str] = set()
    with gzip.open(db, "rt", encoding="latin-1") as fh:
        seq: list[str] = []
        for line in fh:
            if line.startswith(">"):
                _add(("".join(seq)), kmers, k)
                seq = []
            else:
                seq.append(line.strip())
        _add("".join(seq), kmers, k)
    try:
        cache.write_bytes(pickle.dumps(kmers))
    except Exception:
        pass
    return frozenset(kmers)


def _add(seq: str, out: set, k: int) -> None:
    s = seq.upper()
    for i in range(len(s) - k + 1):
        out.add(s[i:i + k])


def predict(peptides: list[str], hosts: list[str], k: int = 9) -> dict[str, dict]:
    """Peptit -> {self_match, host, match}. Herhangi bir konak self-k-mer'iyle tam eşleşme.

    Yumuşak: sadece işaretler. Bir konağın DB'si yoksa o konak atlanır (dürüst).
    """
    host_sets = [(h, _kmers(h, k)) for h in hosts if available(h)]
    host_sets = [(h, s) for h, s in host_sets if s]
    out: dict[str, dict] = {}
    for pep in set(peptides):
        s = pep.upper()
        hit_host, hit_kmer = "", ""
        if len(s) >= k:
            for i in range(len(s) - k + 1):
                sub = s[i:i + k]
                for h, kset in host_sets:
                    if sub in kset:
                        hit_host, hit_kmer = h, sub
                        break
                if hit_kmer:
                    break
        out[pep] = {"self_match": bool(hit_kmer), "host": hit_host, "match": hit_kmer}
    return out
