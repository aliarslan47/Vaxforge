"""Adım 4 (MHC-I sonrası) — antijen işleme: IEDB Analysis Resource 'processing' API.

NetCTL-1.2'NİN YERİNE GEÇEN YOL. netctl.py'nin modül başlığında anlatılan kusur
bu makinede kalıcı: netCTL-1.2a'nın x86_64 binary'si satıcının kendi test
dosyasıyla bile SEGFAULT veriyor (statik bağlı eski binary, yeni çekirdek).
Yerelde düzeltilemiyor → `processing_norm` her adayda None kalıyordu, yani MHC-I
adaylarının sunulabilirliği HİÇ ölçülmüyordu.

IEDB 'processing' aracı aynı iki allel-bağımsız bileşeni HTTP üzerinden verir:
  * proteasome_score — proteozomal C-terminal kesim
  * tap_score        — TAP taşıma yakınlığı

ÖLÇEK UYARISI (önemli): IEDB'nin `proteasome_score`'u NetCTL'in `cle`'si DEĞİLDİR.
NetCTL cle'si 0-1 olasılıktır; IEDB'nin skoru log-ölçekli ve 1'i aşar. Bu yüzden
bu modül KENDİ normalizasyonunu taşır: iki skor da sigmoid ile 0-1'e çekilir,
NetCTL'in kesim:TAP = 3:1 göreli ağırlığı korunur (karşılaştırılabilirlik için).

ÖLÇÜLEN ETKİ (2026-11-17, 105 gerçek IEDB cle değeri, aralık 0.321–1.722):
değerlerin %68'i 1.0'ı aşıyor. netctl.processing_norm() cle'yi [0,1]'e KIRPTIĞI
için bu %68 tek bir değere çöküyor: kesim bileşeni tek başına 105 değerde
35 benzersiz seviyeye düşüyor, sigmoidle 102. Nihai processing_norm yine de
tamamen doymaz (TAP terimi kırpılmadığından değişmeye devam eder) — yani kusur
"her şey 1.0 olur" değil, "kesim çözünürlüğünün üçte ikisi kaybolur". Yorumun
ilk hâli bunu abartıyordu; ölçüp düzeltildi.

DÜRÜSTLÜK: yöntem etiketi "NetCTL-1.2" DEMEZ — farklı tahminci, farklı ölçek.
Ağ yoksa / API hata verirse boş döner (çağıran sessizce atlar, uydurma yok).
Sonuçlar diske önbelleklenir; aynı protein ikinci koşuda tekrar sorulmaz.
"""

from __future__ import annotations

import hashlib
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

_URL = "https://tools-cluster-interface.iedb.org/tools_api/processing/"
_CACHE_DIR = Path(__file__).resolve().parent.parent / "tools" / "db" / "iedb_processing"
_TIMEOUT = 300

METHOD = "GERÇEK (IEDB processing API — proteozomal kesim + TAP; NetCTL-1.2 değil)"
METHOD_EN = "REAL (IEDB processing API — proteasomal cleavage + TAP; not NetCTL-1.2)"

# NetCTL'in içsel göreli ağırlıkları (kesim:TAP = 3:1) — netctl.py ile aynı,
# böylece iki omurganın ürettiği processing_norm aynı ölçekte okunur.
_W_CLE, _W_TAP = 0.15, 0.05

# API allel ister ama proteasome/tap skorları allel-BAĞIMSIZDIR (aynı peptit
# için allel değişince değişmezler); yalnız istek geçerli olsun diye verilir.
# mhc_score/ic50 kolonları OKUNMAZ — MHC bağlanmayı NetMHCpan/IEDB MHC yapıyor.
_DUMMY_ALLELE = "HLA-A*02:01"

# Tek istekte gönderilen protein sayısı. IEDB paylaşımlı bir küme; büyük FASTA
# zaman aşımına ve 502'ye yol açıyor. Küçük parça + önbellek daha güvenilir.
CHUNK = 12
PEPTIDE_LENGTH = 9  # NetCTL de yalnız 9-mer tarar (bkz. netctl.py başlığı)


def available() -> bool:
    """API erişilebilir mi? Önbellek varsa ağ olmadan da True.

    Dürüstlük: ağ yoksa ve önbellek boşsa False — araç 'hazır' iddia etmez.
    """
    if _CACHE_DIR.is_dir() and any(_CACHE_DIR.glob("*.json")):
        return True
    return _probe()


def _probe() -> bool:
    """Tek kısa gerçek çağrı: uç gerçekten hesaplıyor mu?"""
    rows = _post("MKLVGGSAAYDEKLVGGSAAYDEKMNLLRKIVKNRRQ")
    return bool(rows)


def _post(fasta_or_seq: str) -> list[dict]:
    """API'ye POST at, tab-ayrılmış tabloyu satır sözlüklerine çevir."""
    data = urllib.parse.urlencode({
        "method": "netmhcpan",
        "sequence_text": fasta_or_seq,
        "allele": _DUMMY_ALLELE,
        "length": str(PEPTIDE_LENGTH),
        "proteasome": "immuno",
    }).encode()
    try:
        req = urllib.request.Request(_URL, data=data, method="POST")
        with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:
            text = r.read().decode("utf-8", "replace")
    except Exception:
        return []
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or "\t" not in lines[0]:
        return []
    header = lines[0].split("\t")
    if "proteasome_score" not in header or "tap_score" not in header:
        return []  # beklenmeyen şema -> sessizce atla, yanlış kolon okumaktan iyi
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) == len(header):
            rows.append(dict(zip(header, parts)))
    return rows


def processing_norm(cle: float | None, tap: float | None) -> float:
    """IEDB kesim+TAP skorlarını 0-1 işleme skoruna indirger.

    netctl.processing_norm'dan FARKI: cle de sigmoidden geçer, çünkü IEDB'nin
    proteasome_score'u olasılık değil log-ölçekli skordur (bkz. modül başlığı).
    Kırpmak yerine sigmoid kullanmak ayrımı korur.
    """
    if cle is None and tap is None:
        return 0.0
    parts, wsum = 0.0, 0.0
    if cle is not None:
        parts += _W_CLE * (1.0 / (1.0 + math.exp(-cle)))
        wsum += _W_CLE
    if tap is not None:
        parts += _W_TAP * (1.0 / (1.0 + math.exp(-tap)))
        wsum += _W_TAP
    return round(parts / wsum, 4) if wsum else 0.0


def _cache_path(seq: str) -> Path:
    h = hashlib.sha256(seq.encode()).hexdigest()[:32]
    return _CACHE_DIR / f"{h}.json"


def _from_cache(seq: str) -> dict[int, dict[str, float]] | None:
    p = _cache_path(seq)
    if not p.exists():
        return None
    try:
        return {int(k): v for k, v in json.loads(p.read_text()).items()}
    except Exception:
        return None


def _to_cache(seq: str, val: dict[int, dict[str, float]]) -> None:
    try:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        _cache_path(seq).write_text(json.dumps({str(k): v for k, v in val.items()}))
    except Exception:
        pass


def predict(proteins: list[tuple[str, str]]) -> dict[str, dict[int, dict[str, float]]]:
    """[(id, seq)] -> {protein_id: {cterm_pos: {"cle":.., "tap":..}}}.

    netctl.predict ile AYNI dönüş şekli (epitope.py ikisini ayırt etmeden
    okuyabilsin): cterm_pos = 9-mer'in 0-tabanlı C-terminal kalıntı indeksi.
    API 1-tabanlı `end` verir → 0-tabanlıya çevrilir.

    Önbellekte olan protein tekrar sorulmaz; kalanlar CHUNK'lık FASTA'larda
    gönderilir. Bir parça başarısız olursa yalnız o parça eksik kalır.
    """
    if not proteins:
        return {}

    out: dict[str, dict[int, dict[str, float]]] = {}
    pending: list[tuple[str, str]] = []
    for pid, seq in proteins:
        cached = _from_cache(seq)
        if cached is not None:
            out[pid] = cached
        else:
            pending.append((pid, seq))

    for start in range(0, len(pending), CHUNK):
        chunk = pending[start:start + CHUNK]
        # seq_num API'de 1-tabanlı ve FASTA sırasını izler.
        fasta = "".join(f">s{i}\n{seq}\n" for i, (_pid, seq) in enumerate(chunk))
        rows = _post(fasta)
        if not rows:
            continue
        per_seq: dict[int, dict[int, dict[str, float]]] = {}
        for row in rows:
            try:
                sn = int(row["seq_num"])          # 1-tabanlı
                end = int(row["end"])             # 1-tabanlı, kapsayıcı
                cle = float(row["proteasome_score"])
                tap = float(row["tap_score"])
            except (KeyError, ValueError):
                continue
            per_seq.setdefault(sn, {})[end - 1] = {"cle": round(cle, 4),
                                                   "tap": round(tap, 4)}
        for i, (pid, seq) in enumerate(chunk):
            vals = per_seq.get(i + 1)
            if vals:
                out[pid] = vals
                _to_cache(seq, vals)
    return out
