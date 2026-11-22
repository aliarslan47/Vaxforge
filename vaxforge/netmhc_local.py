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

import concurrent.futures
import functools
import os
import platform
import re
import subprocess
import tempfile
from pathlib import Path

from . import mhc_rank

_TOOLS = Path(__file__).resolve().parent.parent / "tools"

# Tek chunk için zaman bütçesi. Maliyet peptit×allel ile büyür; bu makinede
# ölçülen hız ~115 (peptit×allel)/sn → 20.000 iş birimi ≈ 3 dk, TIMEOUT'un
# (20 dk) çok altında. Böylece panel kaç allel olursa olsun chunk timeout'a
# çarpmaz. Üst sınır tek-allelli panelde dosya başına makul kalsın diye var.
_UNITS_PER_CHUNK = 20_000
_MAX_CHUNK = 20_000
TIMEOUT = 1200

# Son predict() çağrısının kapsama tanılaması: {mhc_class: {...}}. Çağıran
# (predictors.py) buradan okuyup yöntem etiketine "kısmi" yazar, böylece
# düşen chunk'lar GERÇEK araç koşmuş gibi görünmez.
_DIAGNOSTIC: dict[str, dict] = {}


def _record(mhc_class: str, **kw) -> None:
    _DIAGNOSTIC[mhc_class] = kw


def last_diagnostic(mhc_class: str) -> dict:
    """Son yerel koşunun kapsaması: requested/returned/failures (boşsa {})."""
    return _DIAGNOSTIC.get(mhc_class, {})

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


def _wrapper_home_ok(wrapper: str | None) -> bool:
    """Sarmalayıcıdaki NMHOME/CTLHOME satırı GERÇEK bir dizini mi gösteriyor?

    DTU paketleri kurulum sırasında elle düzenlenmek üzere satıcı YER TUTUCU
    yoluyla gelir (ör. `setenv NMHOME /tools/src/netMHCpan-4.2`). Düzenlenmezse
    araç hiç çalışmaz ama sessizce başarısız olur. Burada yol bir dizin değilse
    aracı yok sayarız; böylece pipeline dürüst 'vekil' etiketiyle devam eder,
    'GERÇEK araç koştu' yanılgısı oluşmaz.
    """
    if not wrapper:
        return False
    try:
        txt = Path(wrapper).read_text(errors="ignore")
    except Exception:
        return False
    m = re.search(r"^\s*setenv\s+(?:NMHOME|CTLHOME)\s+(\S+)", txt, re.M)
    if not m:
        return True          # bu desen yoksa karar verme (TMHMM/SignalP gibi)
    return Path(m.group(1)).is_dir()


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
    if not _wrapper_home_ok(w):
        return False
    try:
        r = subprocess.run([w, "-h"], capture_output=True, timeout=30, text=True)
        out = (r.stdout + r.stderr).lower()
        # DTU sarmalayıcısı yapılandırılmamışsa "no binaries found for <platform> <yol>"
        # yazar ve hiçbir şey hesaplamaz. Bu metin eskiden kontrol EDİLMİYORDU, bu
        # yüzden runnable() True dönüyor, predict() ise sessizce heuristik vekile
        # düşüyordu: koşu "GERÇEK araç" sanılıp vekil sayılarla tamamlanıyordu.
        for bad in ("format error", "no binaries found", "cannot find", "command not found"):
            if bad in out:
                return False
        return True
    except Exception:
        return False


def available() -> bool:
    return runnable("mhc_i") or runnable("mhc_ii")


@functools.lru_cache(maxsize=2)
def _supported(mhc_class: str) -> frozenset:
    """Aracın desteklediği allel adları (tam araç biçimi, ör. BoLA-DRB3_01101).

    netMHCpan '-listMHC', netMHCIIpan '-list'. Alınamazsa boş küme → çağıran
    doğrulamayı ATLAR (eski davranış, regresyon yok).
    """
    w = _wrapper(mhc_class)
    if not w:
        return frozenset()
    # ONCE ARACIN KENDI VERI DOSYASI: '-listMHC' cikti es-anlamli adlari da basar
    # (SLA-10101 gibi) ama '-a' bayragi YALNIZ data/allelenames (MHC-I) ya da
    # data/allelelist.txt (MHC-II) icerigini kabul eder. Dogrulamayi genis listeye
    # yaparsak gecersiz allel suzgecten gecer, arac hata verir, konak sessizce
    # 0 baglanan alir. Bu yuzden zorlayici liste once okunur.
    names = _enforced_allele_names(w, mhc_class)
    if names:
        return names
    flag = "-listMHC" if mhc_class == "mhc_i" else "-list"
    try:
        r = subprocess.run([w, flag], capture_output=True, timeout=60, text=True)
    except Exception:
        return frozenset()
    names = set()
    for line in r.stdout.splitlines():
        t = line.strip()
        if t and not t.startswith("#"):
            names.add(t.split()[0])
    return frozenset(names)


def _tool_home(wrapper: str | None) -> Path | None:
    """Sarmalayicidaki NMHOME satirindan aracin kurulum dizinini okur."""
    if not wrapper:
        return None
    try:
        txt = Path(wrapper).read_text(errors="ignore")
    except Exception:
        return None
    m = re.search(r"^\s*setenv\s+(?:NMHOME|CTLHOME)\s+(\S+)", txt, re.M)
    if not m:
        return None
    h = Path(m.group(1))
    return h if h.is_dir() else None


@functools.lru_cache(maxsize=2)
def _enforced_allele_names(wrapper: str, mhc_class: str) -> frozenset:
    """'-a' bayraginin GERCEKTEN kabul ettigi allel adlari (zorlayici liste).

    netMHCpan '-a' dogrulamasini data/allelenames, netMHCIIpan ise
    data/allelelist.txt dosyasina karsi yapar. '-listMHC'/'-list' ciktisi bunun
    UST KUMESIDIR (es-anlamli yazimlari da basar), bu yuzden dogrulama icin
    guvenilmez. Dosya bulunamazsa bos kume -> cagiran CLI listesine duser.
    """
    home = _tool_home(wrapper)
    if home is None:
        return frozenset()
    fname = "allelenames" if mhc_class == "mhc_i" else "allelelist.txt"
    for cand in (home / "Linux_x86_64" / "data" / fname, home / "data" / fname):
        if cand.is_file():
            names = set()
            try:
                for line in cand.read_text(errors="ignore").splitlines():
                    t = line.strip()
                    if t and not t.startswith("#"):
                        names.add(t.split()[0])
            except Exception:
                return frozenset()
            return frozenset(names)
    return frozenset()


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
        # HLA tipi adlarda '*' ardindan ZATEN kolon var (HLA-A*02:01) -> '*' atilir
        # ve 'HLA-A02:01' dogru bicimi cikar. Ama SLA/BoLA tipi adlarda kolon YOK
        # (SLA-1*0101); '*' atilinca haneler kaynasir -> 'SLA-10101'. Bu ad araçın
        # '-listMHC' listesinde (es-anlamli olarak) GORUNUR ama gercek kapi olan
        # data/allelenames dosyasinda YOKTUR -> '-a' cagrisi komple duser ve konak
        # SESSIZCE 0 baglanan alir (kuduz kosusunda domuz: 19 -> 0).
        # Dogru bicim 'SLA-1:0101'. Bu yuzden kolon yoksa '*' -> ':' cevrilir.
        star = a.find("*")
        if star >= 0 and ":" not in a[star:]:
            return a.replace("*", ":", 1)
        return a.replace("*", "")
    x = a
    if x.upper().startswith("HLA-DR"):
        x = x[4:]                       # 'HLA-' at -> DRB1*01:01
    return x.replace("*", "_").replace(":", "")


def predict(peptides: list[str], alleles: list[str], mhc_class: str,
            rank_weak: float, rank_strong: float | None = None) -> dict[str, dict] | None:
    """Yerel binary ile tahmin. Çalışmazsa None (çağıran IEDB'ye düşer)."""
    if not runnable(mhc_class):
        return None
    w = _wrapper(mhc_class)
    peps = sorted({p for p in peptides})
    by_len: dict[int, list[str]] = {}
    for p in peps:
        by_len.setdefault(len(p), []).append(p)
    alleles_n = [_norm_allele(a, mhc_class) for a in alleles]
    # Geçersiz allelleri ele: netMHCpan/IIpan -a listesinde TEK geçersiz allel
    # tüm çağrıyı düşürür → hepsi proxy'ye kaçardı. Araç -list'ine karşı süz.
    supported = _supported(mhc_class)
    if supported:
        valid = [a for a in alleles_n if a in supported]
        if not valid:
            return None    # hiç geçerli allel yok → çağıran proxy'ye düşer (dürüst)
        alleles_n = valid
    # Çıktıdaki allel adını (araç biçimi) orijinal konak allel adına (IEDB biçimi,
    # ör. HLA-A*02:01) çevirmek için normalize-anahtar haritası — popülasyon
    # kapsamı gerçek allel adlarını gerektirir.
    key2orig = {_akey(a): a for a in alleles}
    # ÖLÇEK GÜVENLİĞİ (bkz. psortb.py aynı sınıf bug): netMHCpan/netMHCIIpan TEK
    # dev batch'te (tam proteom → yüz binlerce–1.37M peptit) ÇÖKÜYOR/timeout'a çarpıyor
    # → yalnız minik bir kısmı dönüyor (%0), çoğu protein MHC-epitopsuz kalıyor → yüzey
    # antijenleri kayboluyor. Bu yüzden peptitleri CHUNK'lara böl; netMHCpan tek-thread
    # olduğundan chunk'ları ÇEKİRDEKLER ARASINDA PARALEL koştur (16 çekirdek → ~12×).
    # Sonuçlar birleştirilir; bir chunk başarısız olsa diğerleri kurtulur (kısmi > hiç).
    #
    # CHUNK SABİT PEPTİT SAYISI OLAMAZ (2026-11-22, BoHV-1 koşusu): araç maliyeti
    # peptit DEĞİL, peptit×allel ile büyür. Bu makinede ölçüldü: netMHCIIpan tek
    # süreçte ~115 (peptit×allel)/sn. Eski sabit 20.000 peptitlik chunk × 11 BoLA
    # alleli = ~32 dk, yani TIMEOUT'un iki katı → her iki chunk da düştü, `agg` boş
    # kaldı, predict None döndü ve çağıran sessizce heuristik vekile indi. Koşu
    # "MHC-II=proxy" ile bitti, oysa araç kurulu ve çalışıyordu (35.993 15-mer).
    # Bu yüzden chunk İŞ MİKTARINA göre bölünür: panel büyüdükçe chunk küçülür.
    n_panel = max(1, len(alleles_n))
    CHUNK = max(250, min(_MAX_CHUNK, _UNITS_PER_CHUNK // n_panel))
    chunks = [(L, group[k:k + CHUNK])
              for L, group in by_len.items()
              for k in range(0, len(group), CHUNK)]
    workers = min(12, max(1, (os.cpu_count() or 4) - 2), len(chunks))
    # Düşen chunk'lar SESSİZ KALMAZ: hangi uzunlukta kaç peptidin neden kaybolduğu
    # buraya yazılır, çağıran `last_diagnostic()` ile okur ve yöntem etiketine
    # "kısmi" yazar. Eskiden `except Exception: return ""` tek kelime etmeden
    # yutuyordu; bütün batch düştüğünde bile kullanıcı yalnız "proxy" görüyordu.
    failures: list[str] = []

    def _run_chunk(task):
        L, sub = task
        fpath = None
        scratch = None
        try:
            with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
                fh.write(_fasta(sub))
                fpath = fh.name
            cmd = [w, "-f", fpath, "-a", ",".join(alleles_n)]
            if mhc_class == "mhc_i":
                cmd += ["-l", str(L)]
            # EŞZAMANLILIK KİLİDİ FIX'i (2026-10-01, brucella deadlock):
            # 12 paralel netMHCpan örneği ORTAK scratch/tmp'de çakışıp birbirini
            # pipe_read'de kilitliyordu → her takılan chunk timeout'a kadar bir
            # worker'ı bloke ediyordu. Çözüm: her chunk'a İZOLE TMPDIR + makul
            # timeout (takılan chunk hızlı bırakılır; kod kısmi sonucu tolere eder).
            scratch = tempfile.mkdtemp(prefix="netmhc_")
            env = dict(os.environ, TMPDIR=scratch, TMP=scratch, TEMP=scratch)
            r = subprocess.run(cmd, capture_output=True, timeout=TIMEOUT, text=True,
                               env=env, stdin=subprocess.DEVNULL)
            if r.returncode != 0:
                failures.append(f"{L}-mer × {len(sub)} peptit: çıkış kodu "
                                f"{r.returncode} ({(r.stderr or '').strip()[:120]})")
                return ""
            return r.stdout
        except subprocess.TimeoutExpired:
            failures.append(f"{L}-mer × {len(sub)} peptit: {TIMEOUT}s timeout")
            return ""      # bu chunk atlandı; diğerleri devam
        except Exception as exc:
            failures.append(f"{L}-mer × {len(sub)} peptit: {type(exc).__name__}")
            return ""
        finally:
            if fpath:
                try:
                    os.unlink(fpath)
                except OSError:
                    pass
            if scratch:
                import shutil as _sh
                _sh.rmtree(scratch, ignore_errors=True)

    agg: dict[str, dict] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        for stdout in ex.map(_run_chunk, chunks):
            if stdout:
                _parse_into(stdout, agg, alleles_n, rank_weak, mhc_class,
                            rank_strong=rank_strong)
    _record(mhc_class, requested=len(peps), returned=len(agg),
            n_chunks=len(chunks), panel=n_panel, chunk=CHUNK, failures=failures)
    if not agg:
        return None
    for pep, d in agg.items():
        d["coverage_frac"] = d["n_alleles"] / max(1, len(alleles_n))
        d["panel_size"] = len(alleles_n)
        d["bound_alleles"] = sorted({key2orig.get(k, k) for k in d.pop("_bound_keys", set())})
    return agg


def _parse_into(stdout: str, agg: dict, alleles, rank_weak: float, mhc_class: str,
                rank_strong: float | None = None) -> None:
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
        # ÖLÇEK: skor %rank'tan TEK kaynakta türetilir (mhc_rank). Eskiden burada
        # `1 - rank/100` vardı; %rank'ın biyolojik aralığını (0.01-2) 0.98-1.00'e
        # sıkıştırdığı için `mhc_score` pratikte sabitti ve adaylık puanında
        # güçlü/sınırda bağlayıcı ayırt edilemiyordu (ölçüldü: yayılım 0.0196).
        score = mhc_rank.rank_to_score(rank, mhc_class, rank_strong, rank_weak)
        if score > d["score"]:
            d["score"] = round(score, 4)
        if rank < d["rank"]:
            d["rank"] = round(rank, 3)
            d["best_allele"] = allele
        if rank <= rank_weak:
            d["n_alleles"] += 1
            d["_bound_keys"].add(_akey(allele))
