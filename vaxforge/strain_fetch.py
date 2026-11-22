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
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_TIMEOUT = 60
_UA = "VaxForge/1.0 (strain-fetch)"
_CACHE = Path(__file__).resolve().parent.parent / "tools" / "db" / ".strain_auto"

#: Virüs izolatında "tam proteom" sayılmak için en az protein sayısı. Virüs
#: kayıtlarının çoğu tek gen/parça (protein_count 1-20); konservasyon için
#: bunlar işe yaramaz, ortolog bulunamaz ve suş sayısını yalancı şişirirler.
_MIN_VIRUS_PROTEINS = 50
#: Ardışık efetch istekleri arası bekleme (NCBI nezaket sınırı: API anahtarsız
#: saniyede 3 istek).
_EUTILS_DELAY = 0.4


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


def _virus_isolates(taxid: str, limit: int = 40,
                    min_proteins: int = _MIN_VIRUS_PROTEINS) -> list[tuple[str, str]]:
    """Virüs izolatları: (nucleotide accession, izolat adı) — TAM proteomlular.

    NEDEN AYRI YOL (2026-11-22, ölçüldü): `genome/taxon/.../dataset_report`
    (assembly yolu) virüslerde yalnız RefSeq REFERANSINI döndürür — BoHV-1 için
    tam olarak 1 kayıt (`GCF_008777455.1`). Yani "3 suş çek" isteği 1 "suş"
    getiriyordu ve o da girdinin kendisiydi; konservasyon kendisiyle
    karşılaştırmadan %100 çıkıyordu. Virüs izolatları başka uç noktada
    (`virus/taxon/.../dataset_report`): BoHV-1'de 854 kayıt, 71'i tam proteom.

    Parça kayıtlar (tek gen, protein_count 1-20) DIŞLANIR: ortolog taraması
    onlarda boş döner, ama suş sayısını yalancı şişirirler.
    """
    out: list[tuple[str, str]] = []
    token = None
    for _page in range(5):
        params = {"page_size": "200"}
        if token:
            params["page_token"] = token
        url = f"{_API}/virus/taxon/{taxid}/dataset_report?" + urllib.parse.urlencode(params)
        try:
            rep = _get(url)
        except Exception:
            break
        for r in rep.get("reports") or []:
            if not r.get("is_annotated"):
                continue
            if (r.get("protein_count") or 0) < min_proteins:
                continue
            acc = r.get("accession")
            if not acc:
                continue
            name = (r.get("isolate") or {}).get("name") or ""
            out.append((acc, name))
            if len(out) >= limit:
                return out
        token = rep.get("next_page_token")
        if not token:
            break
    return out


def _efetch_protein_faa(accession: str) -> str | None:
    """Nükleotid accession → CDS çevirileri (protein FASTA), E-utilities ile.

    Virüs izolatları için Datasets indirme uç noktası protein FASTA döndürmüyor
    (zip yalnız data_report içeriyor — denendi), bu yüzden efetch
    `rettype=fasta_cds_aa` kullanılıyor.
    """
    params = {"db": "nuccore", "id": accession,
              "rettype": "fasta_cds_aa", "retmode": "text"}
    try:
        txt = _get(f"{_EUTILS}/efetch.fcgi?" + urllib.parse.urlencode(params),
                   binary=True).decode("utf-8", "replace")
    except Exception:
        return None
    return txt if txt.count(">") else None


def _fetch_virus_strains(taxid: str, n: int, exclude: set[str]) -> dict:
    """Virüs yolu: tam proteomlu izolatları efetch ile indir (önbellekli)."""
    isolates = _virus_isolates(taxid, limit=max(20, n * 4))
    if not isolates:
        return {"ok": False, "paths": [], "accessions": [], "isolates": [],
                "note": "bu taxon için tam proteomlu virüs izolatı bulunamadı"}
    paths, accs, names = [], [], []
    for acc, name in isolates:
        if len(paths) >= n:
            break
        if acc in exclude or acc.split(".")[0] in exclude:
            continue
        cache_path = _CACHE / f"{acc}.faa"
        if not (cache_path.exists() and cache_path.stat().st_size > 0):
            faa = _efetch_protein_faa(acc)
            time.sleep(_EUTILS_DELAY)
            if not faa:
                continue
            try:
                cache_path.write_text(faa)
            except Exception:
                continue
        paths.append(str(cache_path))
        accs.append(acc)
        names.append(name or acc)
    if not paths:
        return {"ok": False, "paths": [], "accessions": [], "isolates": [],
                "note": "virüs izolat proteomları indirilemedi"}
    return {"ok": True, "paths": paths, "accessions": accs, "isolates": names,
            "note": f"{len(paths)} virüs izolatı NCBI'den çekildi "
                    f"(tam proteom, >= {_MIN_VIRUS_PROTEINS} protein)"}


def fetch_strains(taxon: str | None, n: int = 3,
                  exclude_accession: str | None = None,
                  profile: str | None = None) -> dict:
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
    exclude = {exclude_accession, (exclude_accession or "").split(".")[0]} - {None, ""}

    # VİRÜS: assembly yolu virüslerde yalnız referansı verir (ölçüldü: BoHV-1'de
    # 1 kayıt) → izolat uç noktasından git. Profil verilmediyse de denenir ve
    # sonuç boşsa assembly yoluna düşülür; yanlış profilde zarar yok.
    if profile == "virus":
        res = _fetch_virus_strains(taxid, n, exclude)
        if res.get("ok"):
            return res
        virus_note = res.get("note")
    else:
        virus_note = None

    assemblies = _list_assemblies(taxid, limit=max(12, n * 3))
    if not assemblies:
        if profile != "virus":
            res = _fetch_virus_strains(taxid, n, exclude)
            if res.get("ok"):
                return res
        return {"ok": False, "paths": [], "accessions": [],
                "note": virus_note or "bu taxon için RefSeq suş bulunamadı"}
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
