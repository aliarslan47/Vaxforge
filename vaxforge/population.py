"""Adım 5b — Popülasyon kapsamı (IEDB Population Coverage, Bui ve ark. 2006).

GERÇEK araç: tools/population-coverage-3.0.2 (IEDB standalone; allelefrequencies.net
insan HLA frekansları paketli). Seçilen aday epitop setinin, farklı coğrafi/etnik
insan popülasyonlarında MHC allel frekanslarına göre nüfusun ne kadarını kapsadığını
hesaplar (bir bireyin en az bir epitop-bağlayan allele sahip olma olasılığı,
Hardy-Weinberg varsayımıyla).

TÜM KONAKLARA GENİŞLETME (kullanıcı kararı): araç her seçili konak için çalıştırılır.
ANCAK IEDB frekans verisi YALNIZ insan HLA'sı içindir — insan gerçek bölgesel kapsam %
verir; insan-dışı konaklarda (fare/sığır/domuz/tavuk) frekans verisi olmadığından araç
sonuç bulamaz ve dürüstçe 'veri yok — hesaplanamadı' döner (uydurma frekans YASAK).
"""

from __future__ import annotations

import functools
import subprocess
import sys
import tempfile
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent.parent / "tools"

METHOD = "GERÇEK (IEDB Population Coverage, Bui 2006)"

# Rapor edilecek temsili coğrafi bölgeler (aracın KENDİ bölge adları — gerçek
# allelefrequencies.net verisi; uydurma değil). World + kıtasal kırılım.
AREAS = ["World", "Europe", "East Asia", "South Asia", "Southeast Asia",
         "North Africa", "West Africa", "North America", "South America", "Oceania"]

_CLASS = {"mhc_i": "I", "mhc_ii": "II"}


def _script() -> Path | None:
    hits = sorted(_TOOLS.glob("population-coverage-*/calculate_population_coverage.py"))
    return hits[0] if hits else None


@functools.lru_cache(maxsize=1)
def available() -> bool:
    return _script() is not None


def _run(pep_alleles: dict[str, list[str]], mhc_class: str,
         areas: list[str]) -> dict[str, dict[str, float]]:
    """{peptit: [allel,...]} -> {bölge: {coverage, average_hit, pc90}}.

    Hiç tanınan allel yoksa (insan-dışı) araç sonuç bulamaz -> boş döner.
    """
    script = _script()
    if script is None or not pep_alleles:
        return {}
    rows = [f"{pep}\t{','.join(al)}" for pep, al in pep_alleles.items() if al]
    if not rows:
        return {}
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write("\n".join(rows) + "\n")
        fpath = fh.name
    try:
        r = subprocess.run(
            [sys.executable, script.name, "-p", ",".join(areas),
             "-c", _CLASS.get(mhc_class, "I"), "-f", fpath],
            capture_output=True, text=True, timeout=300, cwd=str(script.parent))
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return {}
    finally:
        Path(fpath).unlink(missing_ok=True)
    return _parse_summary(r.stdout)


def _parse_summary(stdout: str) -> dict[str, dict[str, float]]:
    """Özet bloğunu ayrıştırır: 'population/area  coverage  average_hit  pc90'
    başlığından sonra bölge satırları ('average' / boş satıra kadar)."""
    out: dict[str, dict[str, float]] = {}
    in_summary = False
    for line in stdout.splitlines():
        s = line.strip()
        if s.startswith("population/area") and "coverage" in s:
            in_summary = True
            continue
        if in_summary:
            if not s or s.startswith(("average", "standard_deviation")):
                break
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            area = parts[0].strip()
            try:
                cov = float(parts[1].replace("%", "").strip())
                avg = float(parts[2].strip())
                pc90 = float(parts[3].strip())
            except ValueError:
                continue
            out[area] = {"coverage": cov, "average_hit": avg, "pc90": pc90}
    return out


def for_candidates(peptides, hosts, areas: list[str] | None = None) -> dict:
    """Aday epitop seti için konak×sınıf×bölge popülasyon kapsamı.

    Her konak için, o konakta bağlanan allellerle {peptit: [allel]} haritası
    kurulur (MHC-I ve MHC-II ayrı) ve IEDB aracıyla koşulur. İnsan gerçek %
    verir; insan-dışı konaklarda araç sonuç bulamayınca 'note' ile dürüstçe
    işaretlenir.
    """
    areas = areas or AREAS
    result = {"method": METHOD, "areas": areas, "available": available(), "hosts": {}}
    if not available():
        result["note"] = "IEDB Population Coverage aracı kurulu değil (tools/)."
        return result

    for h in hosts:
        entry = {"label": h.label}
        any_data = False
        for cls, kind in (("mhc_i", "MHC-I"), ("mhc_ii", "MHC-II")):
            pep_alleles: dict[str, list[str]] = {}
            for p in peptides:
                if p.kind != kind:
                    continue
                per_host = p.metrics.get("per_host", {}) or {}
                bound = (per_host.get(h.name, {}) or {}).get("bound_alleles", [])
                # yalnız insan HLA adları IEDB frekans tablosunda var
                bound = [a for a in bound if a.upper().startswith("HLA")]
                if bound:
                    pep_alleles[p.seq] = bound
            cov = _run(pep_alleles, cls, areas) if pep_alleles else {}
            entry[cls] = cov or None
            any_data = any_data or bool(cov)
        if not any_data:
            entry["note"] = ("IEDB'de bu konak için HLA frekans verisi yok "
                             "(popülasyon kapsamı yalnız insan HLA için hesaplanır).")
        result["hosts"][h.name] = entry
    return result
