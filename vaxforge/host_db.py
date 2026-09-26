"""Konak proteom DB otomatik sağlama — otoimmünite (self-homoloji) için.

funnel'daki ZORUNLU otoimmünite filtresi, seçilen konağın proteom DB'sine karşı
diamond blastp koşar. Kurulu konaklar: insan/fare/sığır/domuz/tavuk. Bunların
dışında bir konak seçilirse DB yoktu → filtre atlanıyordu (korumasız).

Bu modül: DB'si olmayan konak için UniProt'tan reviewed (Swiss-Prot) proteomunu
indirir ve diamond DB kurar (tools/db/{host}.dmnd). Böylece otoimmünite filtresi
her konak için çalışabilir. Ağ/diamond yok → None (dürüst; filtre 'atlandı' der).

Kaynak: UniProt REST (rest.uniprot.org). Swiss-Prot (reviewed) = kompakt, kürasyonlu.
"""

from __future__ import annotations

import shutil
import socket
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent.parent / "tools"
_DB_DIR = _TOOLS / "db"

# Konak adı → UniProt/NCBI taxonomy id (yaygın konaklar). Host.taxid varsa o önceliklidir.
_TAXID = {
    "human": 9606, "mouse": 10090, "bovine": 9913, "cattle": 9913,
    "pig": 9823, "swine": 9823, "chicken": 9031, "rabbit": 9986,
    "rat": 10116, "horse": 9796, "sheep": 9940, "goat": 9925,
    "dog": 9615, "cat": 9685, "monkey": 9544, "macaque": 9544,
    "zebrafish": 7955, "salmon": 8030,
}


def _diamond() -> str | None:
    local = _TOOLS / "diamond"
    if local.exists():
        return str(local)
    return shutil.which("diamond")


def db_path(host_name: str) -> Path:
    return _DB_DIR / f"{host_name}.dmnd"


def available(timeout: float = 3.0) -> bool:
    """UniProt REST + diamond var mı? (otomatik DB sağlanabilir mi)."""
    if _diamond() is None:
        return False
    try:
        with socket.create_connection(("rest.uniprot.org", 443), timeout=timeout):
            return True
    except OSError:
        return False


def _resolve_taxid(host_name: str, taxid: int | None) -> int | None:
    if taxid:
        return int(taxid)
    return _TAXID.get(host_name.lower())


def _download_swissprot(taxid: int, timeout: int = 180) -> str | None:
    """UniProt reviewed (Swiss-Prot) FASTA — organism_id:taxid. Hata → None."""
    q = f"(organism_id:{taxid}) AND (reviewed:true)"
    url = ("https://rest.uniprot.org/uniprotkb/stream?"
           + urllib.parse.urlencode({"format": "fasta", "query": q, "compressed": "false"}))
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "VaxForge/1.0 (host-db)"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", "replace")
    except Exception:
        return None


def ensure(host_name: str, taxid: int | None = None) -> dict:
    """Konak DB'sini sağla. Döner: {ok, path(rel), note, n_proteins, provisioned}.

    Zaten varsa provisioned=False. Yoksa UniProt'tan indirip diamond DB kurar.
    path, funnel'ın beklediği göreli biçimdir (db/{host}.dmnd).
    """
    dmnd = db_path(host_name)
    rel = f"db/{host_name}.dmnd"
    if dmnd.exists():
        return {"ok": True, "path": rel, "note": "mevcut", "provisioned": False}
    dia = _diamond()
    if dia is None:
        return {"ok": False, "path": None, "note": "diamond yok", "provisioned": False}
    tid = _resolve_taxid(host_name, taxid)
    if not tid:
        return {"ok": False, "path": None,
                "note": f"'{host_name}' için taxid bilinmiyor — DB sağlanamadı", "provisioned": False}
    if not available():
        return {"ok": False, "path": None, "note": "UniProt erişilemedi (çevrimdışı)", "provisioned": False}
    faa = _download_swissprot(tid)
    if not faa or faa.count(">") == 0:
        return {"ok": False, "path": None,
                "note": f"UniProt Swiss-Prot boş döndü (taxid {tid})", "provisioned": False}
    _DB_DIR.mkdir(parents=True, exist_ok=True)
    faa_path = _DB_DIR / f"{host_name}_sprot_auto.fasta"
    try:
        faa_path.write_text(faa)
        subprocess.run([dia, "makedb", "--in", str(faa_path), "-d", str(dmnd), "--quiet"],
                       capture_output=True, timeout=300, check=True)
    except Exception as e:
        return {"ok": False, "path": None,
                "note": f"diamond makedb başarısız ({type(e).__name__})", "provisioned": False}
    return {"ok": True, "path": rel, "note": f"UniProt'tan sağlandı (taxid {tid})",
            "n_proteins": faa.count(">"), "provisioned": True}
