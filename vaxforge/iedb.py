"""GERÇEK NetMHCpan / NetMHCIIpan — IEDB Analysis Resource REST API üzerinden.

DTU'nun NetMHCpan/NetMHCIIpan binary'leri akademik lisans ister ve script'le
indirilemez. IEDB API bu algoritmaları ücretsiz, çok-türlü, HTTP üzerinden sunar
(insan HLA + fare H-2 doğrulandı; sığır/tavuk/domuz adlandırması desteklenmiyorsa
predictors.py proxy'ye düşer).

Doğrulanan alleller diske önbelleklenir (tekrar tekrar problamamak için).
İnternet yoksa/hatada boş döner -> çağıran proxy'ye düşer.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

MHCI_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
MHCII_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
METHOD_I = "netmhcpan_el"
METHOD_II = "netmhciipan_el"
_CACHE = Path(__file__).resolve().parent.parent / "outputs" / ".iedb_allele_cache.json"


def _post(url: str, data: dict, timeout: int = 120) -> str | None:
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, data=body, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", "replace")
    except Exception:
        return None


def _parse(text: str | None) -> list[dict] | None:
    if not text or "Invalid" in text[:200] or "error" in text[:80].lower():
        return None
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or "\t" not in lines[0]:
        return None
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) == len(header):
            rows.append(dict(zip(header, parts)))
    return rows or None


# --- allel doğrulama (önbellekli) ------------------------------------------
def _load_cache() -> dict:
    if _CACHE.exists():
        try:
            return json.loads(_CACHE.read_text())
        except Exception:
            return {}
    return {}


def _save_cache(c: dict) -> None:
    _CACHE.parent.mkdir(parents=True, exist_ok=True)
    _CACHE.write_text(json.dumps(c, indent=0))


def valid_alleles(alleles: list[str], mhc_class: str) -> list[str]:
    """Verilen allellerden IEDB'nin desteklediklerini döndürür (önbellekli probe)."""
    cache = _load_cache()
    key = "mhci" if mhc_class == "mhc_i" else "mhcii"
    cache.setdefault(key, {})
    url = MHCI_URL if mhc_class == "mhc_i" else MHCII_URL
    method = METHOD_I if mhc_class == "mhc_i" else METHOD_II
    probe_seq = "SLYNTVATLYCVHQRIDVKDF"
    changed = False
    ok = []
    for a in alleles:
        if a in cache[key]:
            if cache[key][a]:
                ok.append(a)
            continue
        data = {"method": method, "sequence_text": probe_seq, "allele": a}
        if mhc_class == "mhc_i":
            data["length"] = "9"
        valid = _parse(_post(url, data, timeout=60)) is not None
        cache[key][a] = valid
        changed = True
        if valid:
            ok.append(a)
    if changed:
        _save_cache(cache)
    return ok


# --- tahmin ----------------------------------------------------------------
def _fasta(peptides: list[str]) -> str:
    return "".join(f">p{i}\n{p}\n" for i, p in enumerate(peptides))


def predict(peptides: list[str], alleles: list[str], mhc_class: str,
            rank_weak: float, max_peptides: int = 400) -> dict[str, dict]:
    """Peptit -> {score, rank, best_allele, n_alleles, panel_size, coverage_frac}.

    Desteklenen allel yoksa / hata olursa boş döner (çağıran proxy'ye düşer).
    """
    alleles = valid_alleles(alleles, mhc_class)
    if not alleles:
        return {}
    peps = sorted({p for p in peptides})[:max_peptides]
    url = MHCI_URL if mhc_class == "mhc_i" else MHCII_URL
    method = METHOD_I if mhc_class == "mhc_i" else METHOD_II
    rank_col = "percentile_rank" if mhc_class == "mhc_i" else "rank"

    # uzunluğa göre grupla (MHC-I length parametresi eşleşmeli)
    by_len: dict[int, list[str]] = {}
    for p in peps:
        by_len.setdefault(len(p), []).append(p)

    agg: dict[str, dict] = {}
    for L, group in by_len.items():
        data = {"method": method, "sequence_text": _fasta(group),
                "allele": ",".join(alleles)}
        if mhc_class == "mhc_i":
            data["length"] = ",".join([str(L)] * len(alleles))
        rows = _parse(_post(url, data))
        if not rows:
            continue
        for r in rows:
            pep = r.get("peptide", "")
            try:
                rank = float(r.get(rank_col, "100"))
                score = float(r.get("score", "0"))
            except ValueError:
                continue
            d = agg.setdefault(pep, {"score": 0.0, "rank": 99.0, "best_allele": "",
                                     "n_alleles": 0, "panel_size": len(alleles)})
            if score > d["score"]:
                d["score"] = round(score, 4)
            if rank < d["rank"]:
                d["rank"] = round(rank, 3)
                d["best_allele"] = r.get("allele", "")
            if rank <= rank_weak:
                d["n_alleles"] += 1
    for pep, d in agg.items():
        d["coverage_frac"] = d["n_alleles"] / max(1, d["panel_size"])
    return agg
