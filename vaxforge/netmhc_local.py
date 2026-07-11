"""Yerel NetMHCpan / NetMHCIIpan binary backend (opsiyonel, hızlı, çevrimdışı).

Kullanıcı DTU'dan lisanslı binary'yi indirip tools/ altına koyarsa bu modül
onu kullanır. ÖNEMLİ: binary makinenin mimarisiyle uyumlu olmalı (bu makine
x86_64; arm64 build çalışmaz). runnable() bunu bir test çalıştırmasıyla
doğrular; uyumsuz/çalışmazsa None döner ve predictors IEDB API'ye düşer.

Beklenen yerleşim:
  tools/netMHCpan-4.2/netMHCpan            (tcsh sarmalayıcı; NMHOME ayarlı)
  tools/netMHCIIpan-4.3/netMHCIIpan
Ortam değişkenleriyle de geçersiz kılınabilir: NETMHCPAN, NETMHCIIPAN.
"""

from __future__ import annotations

import functools
import os
import platform
import subprocess
import tempfile
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent.parent / "tools"

# ELF e_machine -> mimari adı
_ELF_MACHINE = {0x3E: "x86_64", 0xB7: "aarch64", 0x28: "arm", 0x03: "x86"}
_ARCH_ALIASES = {"x86_64": {"x86_64", "amd64"}, "aarch64": {"aarch64", "arm64"}}


def _wrapper(mhc_class: str) -> str | None:
    env = os.environ.get("NETMHCPAN" if mhc_class == "mhc_i" else "NETMHCIIPAN")
    if env and Path(env).exists():
        return env
    pat = "netMHCpan-*/netMHCpan" if mhc_class == "mhc_i" else "netMHCIIpan-*/netMHCIIpan"
    hits = sorted(_TOOLS.glob(pat))
    return str(hits[0]) if hits else None


def _elf_arch(path: Path) -> str | None:
    try:
        with open(path, "rb") as fh:
            head = fh.read(20)
        if head[:4] != b"\x7fELF":
            return None
        e_machine = int.from_bytes(head[18:20], "little")
        return _ELF_MACHINE.get(e_machine)
    except Exception:
        return None


def _compiled_binary(mhc_class: str) -> Path | None:
    w = _wrapper(mhc_class)
    if not w:
        return None
    root = Path(w).parent
    for d in root.glob("Linux_*/bin/*"):
        if _elf_arch(d):
            return d
    return None


@functools.lru_cache(maxsize=2)
def runnable(mhc_class: str) -> bool:
    """Yerel binary bu MAKİNEDE çalışır mı? Derlenmiş ELF mimarisini kontrol eder.

    arm64 build'i x86_64 makinede çalıştırmaz (sizin durumunuz). Mimari uyarsa
    ayrıca hızlı bir test koşusu dener.
    """
    binp = _compiled_binary(mhc_class)
    if binp is None:
        return False
    arch = _elf_arch(binp)
    machine = platform.machine().lower()
    if arch is None or machine not in _ARCH_ALIASES.get(arch, {arch}):
        return False   # mimari uyuşmazlığı -> çalışmaz
    w = _wrapper(mhc_class)
    try:
        r = subprocess.run([w, "-h"], capture_output=True, timeout=30, text=True)
        out = (r.stdout + r.stderr).lower()
        if "format error" in out:
            return False
        return True
    except Exception:
        return False


def available() -> bool:
    return runnable("mhc_i") or runnable("mhc_ii")


def _fasta(peptides: list[str]) -> str:
    return "".join(f">p{i}\n{p}\n" for i, p in enumerate(peptides))


def _akey(a: str) -> str:
    """Allel adını biçimden bağımsız eşleştirme anahtarına indirger.

    HLA-A*02:01 ve HLA-A02:01 -> 'A0201'; HLA-DRB1*01:01 ve DRB1_0101 -> 'DRB10101'.
    """
    return a.upper().replace("HLA-", "").replace("*", "").replace("_", "").replace(":", "")


def _norm_allele(a: str, mhc_class: str) -> str:
    """Allel adını yerel araç biçimine çevirir.

    MHC-I  (netMHCpan):    HLA-A*02:01 -> HLA-A02:01  ('*' atılır)
    MHC-II (netMHCIIpan):  HLA-DRB1*01:01 -> DRB1_0101 ('HLA-' at, '*'->'_', ':' sil)
    """
    if mhc_class == "mhc_i":
        return a.replace("*", "")
    x = a
    if x.upper().startswith("HLA-DR"):
        x = x[4:]                       # 'HLA-' at -> DRB1*01:01
    return x.replace("*", "_").replace(":", "")


def predict(peptides: list[str], alleles: list[str], mhc_class: str,
            rank_weak: float) -> dict[str, dict] | None:
    """Yerel binary ile tahmin. Çalışmazsa None (çağıran IEDB'ye düşer)."""
    if not runnable(mhc_class):
        return None
    w = _wrapper(mhc_class)
    peps = sorted({p for p in peptides})
    by_len: dict[int, list[str]] = {}
    for p in peps:
        by_len.setdefault(len(p), []).append(p)
    alleles_n = [_norm_allele(a, mhc_class) for a in alleles]
    # Çıktıdaki allel adını (araç biçimi) orijinal konak allel adına (IEDB biçimi,
    # ör. HLA-A*02:01) çevirmek için normalize-anahtar haritası — popülasyon
    # kapsamı gerçek allel adlarını gerektirir.
    key2orig = {_akey(a): a for a in alleles}
    agg: dict[str, dict] = {}
    try:
        for L, group in by_len.items():
            with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
                fh.write(_fasta(group))
                fpath = fh.name
            cmd = [w, "-f", fpath, "-a", ",".join(alleles_n)]
            if mhc_class == "mhc_i":
                cmd += ["-l", str(L)]
            # Global batch (tüm proteinlerin peptitleri tek çağrıda) → büyük olabilir;
            # süre sınırı geniş tutulur (protein-başına eski 600s yerine).
            r = subprocess.run(cmd, capture_output=True, timeout=3600, text=True)
            os.unlink(fpath)
            _parse_into(r.stdout, agg, alleles_n, rank_weak, mhc_class)
    except Exception:
        return None
    if not agg:
        return None
    for pep, d in agg.items():
        d["coverage_frac"] = d["n_alleles"] / max(1, len(alleles_n))
        d["panel_size"] = len(alleles_n)
        d["bound_alleles"] = sorted({key2orig.get(k, k) for k in d.pop("_bound_keys", set())})
    return agg


def _parse_into(stdout: str, agg: dict, alleles, rank_weak: float, mhc_class: str) -> None:
    """netMHCpan/netMHCIIpan tablosunu ayrıştırır (başlık sütun indeksinden).

    Başlık MHC-I'de normal satır, MHC-II'de '#' ile yorumlu; ikisini de destekler.
    Peptit/rank/allel sütun indekslerini başlıktan bulur, veri satırlarını buna göre
    okur (satır sonundaki 'BindLevel <= SB' fazlalıkları sorun olmaz).
    """
    pep_i = rank_i = mhc_i = None
    for line in stdout.splitlines():
        raw = line.rstrip()
        if not raw.strip():
            continue
        toks = raw.lstrip("#").split()
        low = [t.lower() for t in toks]
        if pep_i is None:
            if "peptide" in low and any("rank" in t for t in low):
                pep_i = low.index("peptide")
                rank_i = next(i for i, t in enumerate(low) if "rank" in t)
                mhc_i = low.index("mhc") if "mhc" in low else (low.index("allele") if "allele" in low else 1)
            continue
        st = raw.strip()
        if st.startswith("#") or st.startswith("-"):
            continue
        parts = st.split()
        if len(parts) <= max(pep_i, rank_i):
            continue
        pep = parts[pep_i]
        if len(pep) < 8 or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in pep):
            continue
        try:
            rank = float(parts[rank_i])
        except ValueError:
            continue
        allele = parts[mhc_i] if mhc_i is not None and mhc_i < len(parts) else (alleles[0] if alleles else "?")
        d = agg.setdefault(pep, {"score": 0.0, "rank": 99.0, "best_allele": "",
                                 "n_alleles": 0, "_bound_keys": set()})
        score = max(0.0, 1 - rank / 100)
        if score > d["score"]:
            d["score"] = round(score, 4)
        if rank < d["rank"]:
            d["rank"] = round(rank, 3)
            d["best_allele"] = allele
        if rank <= rank_weak:
            d["n_alleles"] += 1
            d["_bound_keys"].add(_akey(allele))
