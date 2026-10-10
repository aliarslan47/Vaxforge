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


# ---------------------------------------------------------------------------
# PROTEİN DÜZEYİ — FAO/WHO birincil ölçütü (80-aa pencerede ≥%35 kimlik)
# ---------------------------------------------------------------------------
# 2026-10-10 kalibrasyonu (benchmark/calibrate_filters.py): yukarıdaki 6-mer
# kuralı 9-mer peptitlerde iyi çalışıyor (yanlış-pozitif ~%2-3) ama TAM
# PROTEİNLERDE non-alerjen proteinlerin ~%95'ini işaretliyor — uzun dizide 6-mer
# eşleşmesi şansa kalıyor. Bu, 6-mer ölçütünün literatürde bilinen özgüllük
# sorunudur; Codex/FAO-WHO'nun BİRİNCİL ölçütü 80 amino asitlik kayan pencerede
# ≥%35 kimliktir. MEV konstruktu bir proteindir, dolayısıyla orada bu ölçüt
# kullanılır (bkz. mev.py). Yaklaşım: diamond blastp ile hizalama uzunluğu ≥80 ve
# kimlik ≥%35 olan en iyi isabet.

import shutil          # noqa: E402
import subprocess      # noqa: E402
import tempfile        # noqa: E402

_DMND = _DB.parent / "allergens.dmnd"
_WINDOW = 80
_MIN_ID = 35.0


def _diamond_bin() -> str | None:
    local = Path(__file__).resolve().parent.parent / "tools" / "diamond"
    return str(local) if local.exists() else shutil.which("diamond")


def _ensure_db() -> bool:
    """allergens.dmnd yoksa fasta'dan bir kez kurar (idempotent)."""
    if _DMND.exists():
        return True
    dmnd, src = _diamond_bin(), _DB
    if not dmnd or not src.exists():
        return False
    with tempfile.TemporaryDirectory() as td:
        plain = Path(td) / "allergens.faa"
        with gzip.open(src, "rt") as fi:
            plain.write_text(fi.read(), encoding="utf-8")
        try:
            r = subprocess.run([dmnd, "makedb", "--in", str(plain),
                                "-d", str(_DMND.with_suffix(""))],
                               capture_output=True, timeout=900)
        except Exception:
            return False
    return r.returncode == 0 and _DMND.exists()


def protein_available() -> bool:
    return _diamond_bin() is not None and (_DMND.exists() or _DB.exists())


def predict_protein(seq: str) -> dict:
    """Tam protein/MEV konstruktu için alerjenite.

    -> {allergen, identity, align_len, subject, criterion, method}
    allergen=None  : ölçüt uygulanamadı (diamond/DB yok) — DÜRÜST 'bilinmiyor'.
    """
    out = {"allergen": None, "identity": None, "align_len": None, "subject": "",
           "criterion": f"FAO/WHO: ≥%{_MIN_ID:.0f} kimlik / ≥{_WINDOW} aa pencere",
           "method": "uygulanamadı (diamond veya alerjen DB yok)"}
    s = (seq or "").strip().upper()
    if not s or not _ensure_db():
        return out
    dmnd = _diamond_bin()
    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.faa"
        q.write_text(f">query\n{s}\n", encoding="utf-8")
        try:
            r = subprocess.run(
                [dmnd, "blastp", "-q", str(q), "-d", str(_DMND),
                 "-o", str(Path(td) / "hits.tsv"), "--outfmt", "6",
                 "qseqid", "sseqid", "pident", "length", "evalue",
                 "--max-target-seqs", "50", "--quiet"],
                capture_output=True, timeout=600)
            rows = (Path(td) / "hits.tsv").read_text().splitlines() if r.returncode == 0 else []
        except Exception:
            return out
    best = None
    for ln in rows:
        f = ln.split("\t")
        if len(f) < 4:
            continue
        try:
            pid, alen = float(f[2]), int(f[3])
        except ValueError:
            continue
        if alen >= _WINDOW and pid >= _MIN_ID and (best is None or pid > best[0]):
            best = (pid, alen, f[1])
    out["method"] = "GERÇEK (diamond blastp → UniProt KW-0020 alerjen seti)"
    if best:
        out.update(allergen=True, identity=round(best[0], 1),
                   align_len=best[1], subject=best[2])
    else:
        out.update(allergen=False)
    return out
