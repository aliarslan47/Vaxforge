"""HLA → hastalık / patojenite ilişki anotasyonu (salt-yorumlama).

Konak MHC panelindeki allellere literatürde yerleşik hastalık/ilaç/koruma
ilişkilerini bağlar. Skorlamayı ETKİLEMEZ — yalnız rapor/UI için bağlamsal not.
Kaynak: config/hla_disease.yaml (her satır PubMed PMID + DOI ile doğrulanmış).

Kullanım:
    from vaxforge import hla_disease
    hla_disease.annotate("HLA-B*57:01")   -> [{disease, category, note, pmid, doi, citation}, ...]
    hla_disease.annotate_panel(["HLA-B*57:01", ...])  -> {allel: [...], ...} (yalnız eşleşenler)
"""

from __future__ import annotations

import functools
from pathlib import Path

import yaml

_PATH = Path(__file__).resolve().parent.parent / "config" / "hla_disease.yaml"

# Kategori -> UI/rapor için önem tonu (risk/advers = uyarı, koruyucu = olumlu)
CATEGORY_TONE = {"risk": "warn", "adverse": "warn", "protective": "good"}


@functools.lru_cache(maxsize=1)
def _load() -> dict[str, list[dict]]:
    if not _PATH.exists():
        return {}
    with open(_PATH, "r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    return raw.get("associations", {}) or {}


def available() -> bool:
    """Anotasyon verisi yüklü mü (dosya var ve en az bir kayıt)."""
    return bool(_load())


def annotate(allele: str) -> list[dict]:
    """Bir allel için hastalık ilişkilerini döndürür (yoksa boş liste).

    Panel biçimiyle (ör. 'HLA-B*57:01') birebir eşleşir.
    """
    return list(_load().get(allele, []))


def annotate_panel(alleles: list[str]) -> dict[str, list[dict]]:
    """Bir allel listesindeki YALNIZ eşleşen allelleri anotasyonlarıyla döndürür."""
    data = _load()
    out: dict[str, list[dict]] = {}
    for a in alleles:
        hit = data.get(a)
        if hit:
            out[a] = list(hit)
    return out


def summarize(alleles: list[str]) -> dict:
    """Panel için özet: kategori sayıları + eşleşen allel-anotasyon haritası.

    UI/rapor tek çağrıda tüketebilsin diye. n_risk/n_adverse/n_protective + entries.
    """
    entries = annotate_panel(alleles)
    n = {"risk": 0, "adverse": 0, "protective": 0}
    for lst in entries.values():
        for e in lst:
            c = e.get("category")
            if c in n:
                n[c] += 1
    return {
        "available": available(),
        "n_risk": n["risk"],
        "n_adverse": n["adverse"],
        "n_protective": n["protective"],
        "entries": entries,
    }
