"""Epitop-seviyesi moleküler mimikri — konak self-peptidine tam eşleşme.

Bir aday epitop, KONAĞIN kendi proteinindeki bir parçayla birebir aynıysa
'moleküler mimikri' = otoimmünite / merkezi-tolerans (T-hücre negatif seçilim)
riski taşır (Vaccine Design Ch4: otoimmün çapraz-reaktif epitopları çıkar).

KALİBRASYON (2026-10-10) sonrası otoimmünite güvenliğinin ASIL kapısı burasıdır.
Önceden protein-bütünü homolojisi funnel'da sert filtreydi, bu modül ise yumuşak
uyarı katmanı. Ölçüm bunu tersine çevirdi: protein-bütünü kimlik, epitop-düzeyi
riskin vekili değil (S. suis enolazı insana %49.8 benzer ama 9-mer'lerinin yalnız
%1.9'u self ile birebir aynı — benchmark/results/homology_granularity.json).
Şimdi: konak self ile TAM eşleşen epitop ELENİR (n=1000 kontrolde doğrulanmış IEDB
epitoplarında %1.3 yanlış-pozitif), BLOSUM-yakın eşleşme elenmez, skor cezası alır
(TCR dejenerasyonu). Protein düzeyi homoloji yumuşak: işaret + skor cezası.
k=9 (MHC-I çekirdek uzunluğu; paylaşılan 9-mer bir CD8+ epitobu olabilir),
config'ten ayarlanır. Bu k-mer kuralı tam proteine ASLA uygulanmaz: rastgele
non-toksin proteinlerin %66.3'ünü işaretliyor.

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
_AA = "ACDEFGHIKLMNPQRSTVWY"


@functools.lru_cache(maxsize=1)
def _blosum():
    """BLOSUM62 matrisi (lazy). Yoksa None → kısmi mimikri atlanır (dürüst)."""
    try:
        from Bio.Align import substitution_matrices
        return substitution_matrices.load("BLOSUM62")
    except Exception:
        return None

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


def _near_hit(kmer: str, kset: frozenset, min_sub_blosum: float = 1.0):
    """Hamming-1 (tek KONSERVATİF ikame) ile konak self-k-mer'ine kısmi eşleşme.

    kmer'in her pozisyonunda 19 ikame denenir; komşu host set'te VE ikame BLOSUM62≥
    min_sub_blosum (konservatif, ör. I↔V=3) ise kısmi mimikri. Milyonlarca k-mer'e
    brute-force yerine komşu-üretme → hızlı. Döner: (host_kmer, pos, orig, sub, blosum) | None.
    """
    b = _blosum()
    if b is None:
        return None
    for i, orig in enumerate(kmer):
        for aa in _AA:
            if aa == orig:
                continue
            neighbor = kmer[:i] + aa + kmer[i + 1:]
            if neighbor in kset:
                try:
                    sc = float(b[orig, aa])
                except (KeyError, IndexError):
                    sc = -9.0
                if sc >= min_sub_blosum:
                    return (neighbor, i, orig, aa, sc)
    return None


def predict(peptides: list[str], hosts: list[str], k: int = 9) -> dict[str, dict]:
    """Peptit -> {self_match, host, match, near_match, near_host, near_kmer, near_sub, near_blosum}.

    self_match = konak self-k-mer'iyle TAM eşleşme (sert mimikri).
    near_match = tam eşleşme yoksa, tek KONSERVATİF ikameyle KISMİ mimikri (BLOSUM62,
    çapraz-reaktivite riski). İkisi de YUMUŞAK — epitopu elemez, işaretler.
    DB'si olmayan konak atlanır (dürüst).
    """
    host_sets = [(h, _kmers(h, k)) for h in hosts if available(h)]
    host_sets = [(h, s) for h, s in host_sets if s]
    out: dict[str, dict] = {}
    for pep in set(peptides):
        s = pep.upper()
        hit_host, hit_kmer = "", ""
        near = None
        near_host = ""
        if len(s) >= k:
            for i in range(len(s) - k + 1):
                sub = s[i:i + k]
                for h, kset in host_sets:
                    if sub in kset:
                        hit_host, hit_kmer = h, sub
                        break
                if hit_kmer:
                    break
            # tam eşleşme yoksa → kısmi (Hamming-1 konservatif) ara
            if not hit_kmer:
                for i in range(len(s) - k + 1):
                    sub = s[i:i + k]
                    for h, kset in host_sets:
                        nh = _near_hit(sub, kset)
                        if nh:
                            near, near_host = (sub,) + nh, h
                            break
                    if near:
                        break
        rec = {"self_match": bool(hit_kmer), "host": hit_host, "match": hit_kmer,
               "near_match": near is not None, "near_host": near_host,
               "near_kmer": near[0] if near else "", "near_hostkmer": near[1] if near else "",
               "near_sub": (f"{near[3]}{near[2]+1}{near[4]}" if near else ""),
               "near_blosum": (near[5] if near else None)}
        out[pep] = rec
    return out
