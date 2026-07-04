"""MHC bağlanma yordayıcı katmanı — konağa + sınıfa göre backend seçer.

Tek bir birleşik arayüz: predict(peptides, host, mhc_class, rank_weak).
Backend seçimi konağın config'indeki 'predictors' alanına göre yapılır:

  mhcflurry     -> GERÇEK, yalnız insan MHC-I (kurulu, GPU'suz)
  netmhcpan/... -> GERÇEK çok-tür + MHC-II [henüz kurulu değil -> proxy'ye düşer]
  iedb_api      -> IEDB REST [opsiyonel, internet]
  proxy         -> saf-Python vekil (gerçek yok -> dürüst 'heuristik' etiket)

Hiçbir yere kilitlenmez: seçilen backend yoksa proxy'ye düşer ve method
etiketinde hangi gerçek aracın eksik olduğu yazılır.
"""

from __future__ import annotations

from . import iedb, mhc_real, netmhc_local
from .hosts import Host
from .sequtils import KD, sanitize


def _proxy_score(pep: str) -> float:
    s = sanitize(pep)
    if len(s) < 8:
        return 0.0
    anchors = [s[1], s[-1]]
    anchor_h = sum(max(0, KD.get(c, 0)) for c in anchors) / 2
    core_h = sum(KD.get(c, 0) for c in s) / len(s)
    return round(max(0.0, min(1.0, 0.5 * min(1.0, anchor_h / 4) + 0.5 * (0.5 + core_h / 8))), 3)


def _proxy(peptides, alleles, rank_weak, missing_tool):
    out = {}
    n = max(1, len(alleles))
    method = f"heuristik proxy ({missing_tool or 'gerçek yordayıcı'} yok)"
    for pep in set(peptides):
        sc = _proxy_score(pep)
        rank = round(max(0.05, (1 - sc) * 10), 2)
        n_pres = max(1, round(sc * n)) if rank <= rank_weak else 0
        out[pep] = {"score": sc, "rank": rank, "best_allele": alleles[0] if alleles else "?",
                    "n_alleles": n_pres, "panel_size": n, "coverage_frac": n_pres / n,
                    "method": method}
    return out


def _from_mhcflurry(peptides, alleles, rank_weak, host, mhc_class):
    raw = mhc_real.predict(peptides, alleles=alleles, weak_percentile=rank_weak)
    if not raw:
        return None
    out = {}
    for pep, d in raw.items():
        psize = d.get("panel_size", len(alleles)) or 1
        out[pep] = {"score": d["presentation_score"], "rank": d["best_percentile"],
                    "best_allele": d["best_allele"], "n_alleles": d["n_alleles"],
                    "panel_size": psize, "coverage_frac": d["n_alleles"] / psize,
                    "method": f"GERÇEK (MHCflurry, {host.label} MHC-I çevrimdışı)"}
    return out


def _from_iedb(peptides, alleles, mhc_class, rank_weak, host):
    raw = iedb.predict(peptides, alleles, mhc_class, rank_weak)
    if not raw:
        return None
    tool = "NetMHCpan" if mhc_class == "mhc_i" else "NetMHCIIpan"
    for pep, d in raw.items():
        d["method"] = f"GERÇEK ({tool}, IEDB API, {host.label} {mhc_class.upper()})"
    return raw


def predict(peptides: list[str], host: Host, mhc_class: str, rank_weak: float) -> dict[str, dict]:
    """Peptit -> {score, rank, best_allele, n_alleles, coverage_frac, method}.

    Öncelik: seçilen gerçek backend (NetMHCpan/NetMHCIIpan via IEDB) -> insan MHC-I
    için MHCflurry çevrimdışı yedek -> proxy. Hiçbir yere kilitli değil; kullanılan
    yöntem 'method' etiketinde yazılır.
    """
    alleles = host.alleles(mhc_class)
    backend = host.predictor(mhc_class)
    human_hla = any(a.startswith("HLA") for a in alleles)

    # 1) GERÇEK: NetMHCpan / NetMHCIIpan
    if backend in ("netmhcpan", "netmhciipan"):
        # 1a) yerel binary (varsa ve bu mimaride çalışıyorsa) — en hızlı, çevrimdışı
        local = netmhc_local.predict(peptides, alleles, mhc_class, rank_weak)
        if local:
            tool = "NetMHCpan" if mhc_class == "mhc_i" else "NetMHCIIpan"
            for d in local.values():
                d["method"] = f"GERÇEK ({tool} yerel, {host.label} {mhc_class.upper()})"
            return local
        # 1b) IEDB API (çok-tür, gerçek)
        res = _from_iedb(peptides, alleles, mhc_class, rank_weak, host)
        if res is not None:
            return res
        # 1c) insan MHC-I ise MHCflurry çevrimdışı yedek
        if mhc_class == "mhc_i" and human_hla and mhc_real.available():
            res = _from_mhcflurry(peptides, alleles, rank_weak, host, mhc_class)
            if res is not None:
                return res
        return _proxy(peptides, alleles, rank_weak, backend)

    # 2) MHCflurry doğrudan seçilmişse (insan MHC-I)
    if backend == "mhcflurry" and mhc_class == "mhc_i" and mhc_real.available():
        res = _from_mhcflurry(peptides, alleles, rank_weak, host, mhc_class)
        if res is not None:
            return res
        return _proxy(peptides, alleles, rank_weak, "MHCflurry-desteği")

    return _proxy(peptides, alleles, rank_weak, backend)
