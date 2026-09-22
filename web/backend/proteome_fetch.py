"""Taxon ID'den referans proteomu NCBI Datasets'ten otomatik çeker.

Kullanıcı dosya yüklemek yerine organizmayı seçip "NCBI'den getir" derse:
o organizmanın RefSeq referans assembly'sinin protein FASTA'sı (protein.faa)
indirilir ve pipeline girdisi olarak kullanılır. Redundant değildir — tek
referans proteom (örn. Monkeypox 10244 → GCF_000857045.1, 180 protein).

Kaynak: NCBI Datasets v2 REST API (O'Leary ve ark., NCBI). Sadece stdlib
urllib + zipfile. İndirilen proteom yerelde cache'lenir (tekrar hızlı).

Dönüş sözleşmesi (fetch_reference_proteome):
  {"ok": bool, "path": str|None, "accession": str|None, "organism": str|None,
   "n_proteins": int, "note": str|None}
"""
from __future__ import annotations

import io
import json
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
_TIMEOUT = 60
_UA = "VaxForge/1.0 (proteome-fetch)"
_CACHE = Path(__file__).resolve().parent.parent.parent / "tools" / "db" / "proteomes"


def _norm_taxon(taxon: str) -> str | None:
    """"NCBITaxon:10244" / "10244" → "10244". Geçersizse None."""
    if not taxon:
        return None
    t = taxon.split(":")[-1].strip()
    return t if t.isdigit() else None


def _get(url: str, binary: bool = False):
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:
        data = r.read()
    return data if binary else json.loads(data.decode("utf-8", "replace"))


def _pick_assembly(taxid: str) -> tuple[str | None, str | None]:
    """Taxon için RefSeq assembly accession'ı + organizma adı seçer.

    Öncelik: 'reference genome' kategorisi > ilk RefSeq assembly.
    """
    params = {"filters.assembly_source": "refseq", "page_size": "20"}
    url = f"{_API}/genome/taxon/{taxid}/dataset_report?" + urllib.parse.urlencode(params)
    rep = _get(url)
    reports = rep.get("reports") or []
    if not reports:
        return None, None
    ref = None
    for r in reports:
        cat = (r.get("assembly_info") or {}).get("refseq_category")
        if cat and "reference" in str(cat).lower():
            ref = r
            break
    r = ref or reports[0]
    acc = r.get("accession")
    org = (r.get("organism") or {}).get("organism_name")
    return acc, org


def _download_protein_faa(accession: str) -> str:
    """Assembly'nin protein.faa içeriğini (metin) döndürür."""
    url = (f"{_API}/genome/accession/{accession}/download"
           f"?include_annotation_type=PROT_FASTA")
    blob = _get(url, binary=True)
    z = zipfile.ZipFile(io.BytesIO(blob))
    faa = [n for n in z.namelist() if n.endswith("protein.faa")]
    if not faa:
        raise ValueError("indirilen pakette protein.faa yok (proteom anotasyonu olmayabilir)")
    return z.read(faa[0]).decode("utf-8", "replace")


def fetch_reference_proteome(taxon: str) -> dict:
    """Taxon → referans proteom .faa (yerel dosya yolu). Cache'li."""
    taxid = _norm_taxon(taxon)
    if not taxid:
        return {"ok": False, "path": None, "accession": None, "organism": None,
                "n_proteins": 0, "note": "geçersiz taxon"}

    _CACHE.mkdir(parents=True, exist_ok=True)
    cache_path = _CACHE / f"NCBITaxon_{taxid}.faa"
    meta_path = _CACHE / f"NCBITaxon_{taxid}.json"

    # Cache isabeti
    if cache_path.exists() and meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
            return {"ok": True, "path": str(cache_path),
                    "accession": meta.get("accession"), "organism": meta.get("organism"),
                    "n_proteins": meta.get("n_proteins", 0), "note": "cache"}
        except Exception:
            pass

    try:
        acc, org = _pick_assembly(taxid)
        if not acc:
            return {"ok": False, "path": None, "accession": None, "organism": org,
                    "n_proteins": 0,
                    "note": "bu taxon için RefSeq referans proteomu yok — dosya yükleyin"}
        faa = _download_protein_faa(acc)
        n = faa.count(">")
        if n == 0:
            return {"ok": False, "path": None, "accession": acc, "organism": org,
                    "n_proteins": 0, "note": "proteom boş döndü"}
        cache_path.write_text(faa)
        meta = {"accession": acc, "organism": org, "n_proteins": n, "taxid": taxid}
        meta_path.write_text(json.dumps(meta, ensure_ascii=False))
        return {"ok": True, "path": str(cache_path), "accession": acc,
                "organism": org, "n_proteins": n, "note": None}
    except Exception as exc:
        return {"ok": False, "path": None, "accession": None, "organism": None,
                "n_proteins": 0, "note": f"NCBI'den çekilemedi: {type(exc).__name__}"}
