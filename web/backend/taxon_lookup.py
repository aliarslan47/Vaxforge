"""Organizma adından NCBI Taxonomy ID araması (Entrez E-utilities).

Kullanıcı taxon ID'sini ezbere bilmek zorunda kalmasın diye: bir organizma
adı (örn. "monkeypox") girilir, NCBI'nin resmi E-utilities servisinde
aranır ve eşleşen taxon'lar (id + bilimsel ad + rütbe) döndürülür. Sadece
stdlib urllib kullanılır (proje geneli tercih). İnternet gerekir; yoksa
dürüstçe boş sonuç + hata notu döner.

Kaynak: NCBI Entrez Programming Utilities (Sayers ve ark., NAR 2022).
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request

_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_TIMEOUT = 12
_UA = "VaxForge/1.0 (taxon-lookup)"
# NCBI, tool + email verilmesini önerir (rate-limit nezaketi).
_TOOL = "vaxforge"
_EMAIL = "aliarslan47@gmail.com"


def _get_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:
        return json.loads(r.read().decode("utf-8", "replace"))


def search(query: str, retmax: int = 8) -> dict:
    """Organizma adını NCBI Taxonomy'de arar.

    Dönüş: {"query", "results": [{"taxon_id", "iri", "name", "rank", "division"}],
            "available": bool, "note": str|None}
    """
    q = (query or "").strip()
    if not q:
        return {"query": q, "results": [], "available": True, "note": None}

    # Kullanıcı doğrudan bir ID yazmış olabilir ("10244" veya "NCBITaxon:10244").
    # O zaman ad araması yerine ID'yi esummary ile çözüp adı getir.
    id_candidate = q.split(":")[-1].strip() if q.lower().startswith("ncbitaxon:") else q
    direct_id = id_candidate if id_candidate.isdigit() else None

    common = {"tool": _TOOL, "email": _EMAIL, "retmode": "json"}
    try:
        if direct_id:
            ids = [direct_id]
        else:
            # 1) esearch → taxon UID listesi
            es_params = {**common, "db": "taxonomy", "term": q, "retmax": str(retmax)}
            es = _get_json(f"{_BASE}/esearch.fcgi?" + urllib.parse.urlencode(es_params))
            ids = (es.get("esearchresult") or {}).get("idlist") or []
        if not ids:
            return {"query": q, "results": [], "available": True,
                    "note": "eşleşme yok"}

        # 2) esummary → her UID için bilimsel ad + rütbe
        su_params = {**common, "db": "taxonomy", "id": ",".join(ids)}
        su = _get_json(f"{_BASE}/esummary.fcgi?" + urllib.parse.urlencode(su_params))
        res = su.get("result") or {}
        out = []
        for uid in res.get("uids", ids):
            rec = res.get(uid) or {}
            name = rec.get("scientificname") or rec.get("commonname") or ""
            if not name:  # geçersiz ID → esummary boş/hatalı kayıt döner
                continue
            out.append({
                "taxon_id": uid,
                "iri": f"NCBITaxon:{uid}",
                "name": name,
                "rank": rec.get("rank") or "",
                "division": rec.get("division") or "",
            })
        if not out:
            return {"query": q, "results": [], "available": True, "note": "eşleşme yok"}
        return {"query": q, "results": out, "available": True, "note": None}
    except Exception as exc:  # ağ yok / NCBI erişilemiyor → dürüst hata
        return {"query": q, "results": [], "available": False,
                "note": f"NCBI'ye erişilemedi: {type(exc).__name__}"}
