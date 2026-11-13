"""Adım 4 (MHC-I sonrası) — antijen işleme: proteozomal C-terminal kesim + TAP.

GERÇEK araç: NetCTL-1.2 (Larsen ve ark. 2005, Eur J Immunol 35:2295). Bir
peptidin sunulabilmesi için yalnız MHC-I'e bağlanması yetmez; önce proteozom
tarafından doğru C-terminalden kesilmeli ve TAP ile ER'ye taşınmalıdır. NetCTL
bunun için iki allel-BAĞIMSIZ skor üretir:

  * cle — C-terminal proteozomal kesim olasılığı (0-1, sinir ağı).
  * tap — TAP taşıma yakınlığı (log-odds, ~ -3..+3).

TASARIM KARARI: NetCTL'in KENDİ eski MHC bağlanma tahmincisini (yalnız 12
supertype) KULLANMIYORUZ — MHC bağlanmayı NetMHCpan-4.2 yapıyor. Buradan yalnız
allel-bağımsız `cle` + `tap` alınır ve funnel-survivor proteinlerinin her 9-mer'i
için C-terminal pozisyona göre indekslenir; epitope._mhc bunları MHC-I
peptitlerimize (peptidin C-terminaline denk gelen 9-mer) eşler.

NetCTL yalnız 9-mer tarar. 8/10/11-mer peptitlerimiz için, aynı C-terminale sahip
9-mer'in cle/tap'ı kullanılır: kesim genuinely bir C-terminal sinyalidir (peptit
uzunluğundan bağımsız geçerli), TAP için ise yaklaşıklıktır (dürüstçe böyle
etiketlenir).
"""

from __future__ import annotations

import functools
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent.parent / "tools"

METHOD = "GERÇEK (NetCTL-1.2 kesim+TAP, Larsen 2005)"

# NetCTL'in içsel göreli ağırlıkları: kesim (wc=0.15) : TAP (wt=0.05) = 3:1.
# MHC bileşenini (we=1.0) DIŞARIDA bırakırız (onu NetMHCpan yapıyor).
_W_CLE, _W_TAP = 0.15, 0.05

_ROW = re.compile(
    r"\bID\s+(?P<id>\S+)\s+pep\s+(?P<pep>\S+).*?\bcle\s+(?P<cle>-?\d+\.\d+)"
    r"\s+tap\s+(?P<tap>-?\d+\.\d+)"
)


def _wrapper() -> str | None:
    hits = sorted(_TOOLS.glob("netCTL-*/netCTL")) + sorted(_TOOLS.glob("netctl-*/netCTL"))
    return str(hits[0]) if hits else None


def _home_ok(wrapper: str) -> bool:
    """Sarmalayıcıdaki `setenv CTLHOME <yol>` gerçek bir dizin mi?

    DTU paketleri satıcı YER TUTUCU yoluyla gelir; düzenlenmezse araç hiç
    çalışmaz (bkz. netmhc_local._wrapper_home_ok, aynı kusur sınıfı).
    """
    try:
        txt = Path(wrapper).read_text(errors="ignore")
    except Exception:
        return False
    m = re.search(r"^\s*setenv\s+CTLHOME\s+(\S+)", txt, re.M)
    return True if not m else Path(m.group(1)).is_dir()


@functools.lru_cache(maxsize=1)
def available() -> bool:
    """NetCTL bu ortamda GERÇEKTEN koşabiliyor mu? (duman testiyle)

    NEDEN DUMAN TESTİ: eski kontrol yalnız sarmalayıcı dosyası + tcsh varlığına
    bakıyordu ve True döndürüyordu. Ama netCTL-1.2a'nın x86_64 binary'si bu
    çekirdekte SEGFAULT veriyor (satıcının kendi test.fsa dosyasıyla bile;
    statik bağlı eski binary). Sonuç: `predict()` boş sözlük dönüyor, her
    adayın `processing_norm` değeri None kalıyor, AMA rapor Yöntemler
    tablosunda NetCTL'i atıfla birlikte KULLANILMIŞ gibi listeliyordu.
    Yani "araç yok" durumu, ölçülmemiş bir metrik + yanlış bir yöntem iddiası
    hâline geliyordu. Duman testi bunu dürüst bir yokluğa çevirir.
    """
    w = _wrapper()
    if not w or shutil.which("tcsh") is None or not _home_ok(w):
        return False
    # Minik gerçek çağrı: çökme/boş çıktı -> araç yok say.
    fpath = None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
            fh.write(">smoke\nMKLVGGSAAYDEKLVGGSAAYDEKMNLLRKIVKNRRQ\n")
            fpath = fh.name
        r = subprocess.run([w, fpath], capture_output=True, text=True, timeout=120)
    except Exception:
        return False
    finally:
        if fpath:
            Path(fpath).unlink(missing_ok=True)
    if r.returncode != 0:
        return False
    return any(_ROW.search(ln) for ln in r.stdout.splitlines())


def processing_norm(cle: float | None, tap: float | None) -> float:
    """cle (0-1) + tap (log-odds) -> 0-1 antijen işleme skoru.

    TAP log-odds sigmoid ile 0-1'e çekilir; NetCTL'in kesim:TAP=3:1 göreli
    ağırlığıyla birleştirilir. İkisi de yoksa 0.0.
    """
    if cle is None and tap is None:
        return 0.0
    parts, wsum = 0.0, 0.0
    if cle is not None:
        parts += _W_CLE * max(0.0, min(1.0, cle))
        wsum += _W_CLE
    if tap is not None:
        parts += _W_TAP * (1.0 / (1.0 + math.exp(-tap)))
        wsum += _W_TAP
    return round(parts / wsum, 4) if wsum else 0.0


def predict(proteins: list[tuple[str, str]]) -> dict[str, dict[int, dict[str, float]]]:
    """[(id, seq)] -> {protein_id: {cterm_pos: {"cle":.., "tap":..}}}.

    cterm_pos = 9-mer'in 0-tabanlı C-terminal kalıntı indeksi (start+8). Araç
    yoksa / hata olursa boş döner (çağıran yer sessizce atlar).
    """
    if not available() or not proteins:
        return {}
    w = _wrapper()
    # Sentetik kısa ID (NetCTL uzun/özel ID'leri kısaltabilir) -> gerçek ID eşlemesi.
    id_map = {f"s{i}": pid for i, (pid, _seq) in enumerate(proteins)}
    with tempfile.NamedTemporaryFile("w", suffix=".fsa", delete=False) as fh:
        for i, (_pid, seq) in enumerate(proteins):
            fh.write(f">s{i}\n{seq}\n")
        fpath = fh.name
    try:
        r = subprocess.run([w, fpath], capture_output=True, text=True, timeout=900)
    except Exception:
        Path(fpath).unlink(missing_ok=True)
        return {}
    finally:
        Path(fpath).unlink(missing_ok=True)

    out: dict[str, dict[int, dict[str, float]]] = {}
    # Her protein için satırlar sırayla gelir; pozisyon = o protein içindeki
    # 9-mer sırası (0-tabanlı start). C-terminal = start + 8.
    order: dict[str, int] = {}
    for line in r.stdout.splitlines():
        m = _ROW.search(line)
        if not m:
            continue
        sid = m.group("id")
        pid = id_map.get(sid, sid)
        start = order.get(sid, 0)
        order[sid] = start + 1
        cterm = start + 8
        out.setdefault(pid, {})[cterm] = {
            "cle": round(float(m.group("cle")), 4),
            "tap": round(float(m.group("tap")), 4),
        }
    return out
