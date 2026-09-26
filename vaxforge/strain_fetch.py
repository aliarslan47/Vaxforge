"""Otomatik çok-suş çekme — NCBI Datasets'ten aynı türün N suş proteomu.

Konservasyon analizi için kullanıcı ek suş vermese bile, organizma taxon'u
verilmişse aynı türün birkaç RefSeq assembly'si (suşu) otomatik indirilir.
Böylece epitop-seviyesi korunmuşluk (suşlar-arası) varsayılan olarak hesaplanabilir.

Kaynak: NCBI Datasets v2 REST API (stdlib urllib + zipfile). İndirilen suşlar
tools/db/.strain_auto/ altında önbelleklenir. Ağ yok / assembly yok → boş liste
(dürüst) → conservation 'hesaplanmadı' der (uydurma yok).
"""

from __future__ import annotations

import io
import json
import socket
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
_TIMEOUT = 60
_UA = "VaxForge/1.0 (strain-fetch)"
_CACHE = Path(__file__).resolve().parent.parent / "tools" / "db" / ".strain_auto"


def _norm_taxon(taxon: str | None) -> str | None:
    if not taxon:
        return None
    t = str(taxon).split(":")[-1].strip()
    return t if t.isdigit() else None


def available(taxon: str | None = None, timeout: float = 3.0) -> bool:
    """NCBI Datasets API erişilebilir mi? (taxon verilmişse geçerli olmalı)."""
    if taxon is not None and _norm_taxon(taxon) is None:
        return False
    try:
        with socket.create_connection(("api.ncbi.nlm.nih.gov", 443), timeout=timeout):
            return True
    except OSError:
        return False


def _get(url: str, binary: bool = False):
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:
        data = r.read()
    return data if binary else json.loads(data.decode("utf-8", "replace"))


def _list_assemblies(taxid: str, limit: int = 12) -> list[tuple[str, str]]:
    """Taxon için (accession, organism) listesi — referans önce. Hata → []."""
    params = {"filters.assembly_source": "refseq", "page_size": str(limit)}
    url = f"{_API}/genome/taxon/{taxid}/dataset_report?" + urllib.parse.urlencode(params)
    try:
        rep = _get(url)
    except Exception:
        return []
    out = []
    for r in rep.get("reports") or []:
        acc = r.get("accession")
        org = (r.get("organism") or {}).get("organism_name") or ""
        if acc:
            out.append((acc, org))
    return out


def _download_protein_faa(accession: str) -> str | None:
    try:
        url = (f"{_API}/genome/accession/{accession}/download"
               f"?include_annotation_type=PROT_FASTA")
        blob = _get(url, binary=True)
        z = zipfile.ZipFile(io.BytesIO(blob))
        faa = [n for n in z.namelist() if n.endswith("protein.faa")]
        if not faa:
            return None
        return z.read(faa[0]).decode("utf-8", "replace")
    except Exception:
        return None


def fetch_strains(taxon: str | None, n: int = 3,
                  exclude_accession: str | None = None) -> dict:
    """Taxon → en fazla n suş proteomu (yerel .faa yolları). Cache'li.

    Döner: {"ok": bool, "paths": [str], "accessions": [str], "note": str}.
    Hata/ağ-yok → ok=False, paths=[] (conservation dürüstçe 'hesaplanmadı' der).
    """
    taxid = _norm_taxon(taxon)
    if not taxid:
        return {"ok": False, "paths": [], "accessions": [], "note": "geçersiz/eksik taxon"}
    if not available(taxon):
        return {"ok": False, "paths": [], "accessions": [], "note": "NCBI erişilemedi (çevrimdışı)"}
    _CACHE.mkdir(parents=True, exist_ok=True)
    assemblies = _list_assemblies(taxid, limit=max(12, n * 3))
    if not assemblies:
        return {"ok": False, "paths": [], "accessions": [],
                "note": "bu taxon için RefSeq suş bulunamadı"}
    paths, accs = [], []
    for acc, _org in assemblies:
        if len(paths) >= n:
            break
        if exclude_accession and acc == exclude_accession:
            continue
        cache_path = _CACHE / f"{acc}.faa"
        if cache_path.exists() and cache_path.stat().st_size > 0:
            paths.append(str(cache_path)); accs.append(acc); continue
        faa = _download_protein_faa(acc)
        if not faa or faa.count(">") == 0:
            continue
        try:
            cache_path.write_text(faa)
        except Exception:
            continue
        paths.append(str(cache_path)); accs.append(acc)
    if not paths:
        return {"ok": False, "paths": [], "accessions": [],
                "note": "suş proteomu indirilemedi (anotasyon yok olabilir)"}
    return {"ok": True, "paths": paths, "accessions": accs,
            "note": f"{len(paths)} suş NCBI'den otomatik çekildi"}
